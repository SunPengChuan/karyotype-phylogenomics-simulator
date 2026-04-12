
import os
import re
import copy
import itertools
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from ete3 import Tree, TreeStyle, TextFace

# -----------------------------
# 1. CONFIGURATION
# -----------------------------
CONFIG = dict(
    input_tree="output_simulator/tree.nwk",
    input_karyotypes="output_simulator/karyotypes_fission.txt",
    input_true_root_karyotype="output_simulator/karyotypes_true_root.txt",
    output_dir="output_reconstruction_fission",
    tree_viz_output_dir="output_simulator",
    enable_tree_viz=True,
    enable_root_dotplot=True,
    min_block_size=2,
    enable_residual_rct_discovery=False,
    root_final_pending_min_count=2,
    outgroup_name="Outgroup",
)

# -----------------------------
# 2. DATA LOADER
# -----------------------------
class DataLoader:
    def __init__(self, config):
        self.cfg = config
        self.species_karyotypes = {}
        self.gene_position_map = {}
        self.outgroup_karyotypes = {}

    def load_data(self):
        if not os.path.exists(self.cfg["input_tree"]):
            raise FileNotFoundError(f"Tree file not found: {self.cfg['input_tree']}")
        tree = Tree(self.cfg["input_tree"], format=1)
        if tree.get_tree_root().name != "Root":
            roots = tree.search_nodes(name="Root")
            if roots:
                tree = roots[0].copy(method="deepcopy")
        for leaf in list(tree.iter_leaves()):
            if leaf.name.startswith("Sp"):
                continue
            leaf.detach()
        
        if not os.path.exists(self.cfg["input_karyotypes"]):
            raise FileNotFoundError(f"Karyotype file not found: {self.cfg['input_karyotypes']}")

        df = pd.read_csv(
            self.cfg["input_karyotypes"],
            sep="\t",
            header=None,
            engine="python",
            comment="-",
            dtype=str,
            names=["raw"],
        ).fillna("")

        self.gene_position_map = {}
        self.outgroup_karyotypes = {}
        outgroup_name = self.cfg.get("outgroup_name", "Outgroup")

        for raw in df["raw"].tolist():
            if not raw:
                continue
            if raw.startswith("Species"):
                continue
            parts = str(raw).split()
            if len(parts) < 3:
                continue

            sp_name = parts[0]
            chr_id = parts[1]
            start_idx = 3 if (len(parts) >= 4 and parts[2].isdigit()) else 2
            genes = [g.replace("-", "") for g in parts[start_idx:]]
            
            if sp_name == outgroup_name:
                self.outgroup_karyotypes[chr_id] = genes
                continue
                
            if not sp_name.startswith("Sp"):
                continue
            
            if sp_name not in self.species_karyotypes:
                self.species_karyotypes[sp_name] = {}
                self.gene_position_map[sp_name] = {}
            
            final_id = chr_id
            self.species_karyotypes[sp_name][final_id] = genes
            
            for idx, g in enumerate(genes):
                self.gene_position_map[sp_name][g] = (final_id, idx + 1)

        missing = [n.name for n in tree.iter_leaves() if n.name not in self.species_karyotypes]
        if missing:
            raise ValueError(f"Karyotypes missing leaf species: {', '.join(sorted(missing))}")
            
        return tree, self.species_karyotypes, self.gene_position_map

def _load_true_root_karyotype(path):
    if not path or not os.path.exists(path):
        return None

    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        engine="python",
        comment="-",
        dtype=str,
        names=["raw"],
    ).fillna("")

    karyo = {}
    for raw in df["raw"].tolist():
        if not raw:
            continue
        if raw.startswith("Species"):
            continue
        parts = str(raw).split()
        if len(parts) < 3:
            continue
        sp_name = parts[0]
        if sp_name != "TrueRoot":
            continue
        cid = parts[1]
        start_idx = 3 if (len(parts) >= 4 and parts[2].isdigit()) else 2
        genes = [g.replace("-", "") for g in parts[start_idx:]]
        karyo[cid] = genes
    
    return karyo

def generate_root_vs_trueroot_dotplot(reconstructed_tsv_path, true_root_karyo_path, out_path, root_label="Root"):
    if not os.path.exists(reconstructed_tsv_path):
        return False
    true_root = _load_true_root_karyotype(true_root_karyo_path)
    if not true_root:
        return False

    df = pd.read_csv(reconstructed_tsv_path, sep="\t", dtype=str).fillna("")
    if "Node" not in df.columns or "ChrID" not in df.columns or "Genes" not in df.columns:
        return False
    df = df[df["Node"] == root_label]

    rec_karyo = {}
    rec_order = []
    for row in df.itertuples(index=False):
        cid = getattr(row, "ChrID", "")
        genes_str = getattr(row, "Genes", "")
        if not cid:
            continue
        genes = [g for g in str(genes_str).split(" ") if g]
        rec_karyo[cid] = genes
        rec_order.append(cid)

    if not rec_karyo:
        return False
    def _flatten(karyo, order):
        genes = []
        for cid in order:
            genes.extend(karyo.get(cid, []))
        return genes

    true_order = sorted(true_root.keys(), key=str)
    true_flat = _flatten(true_root, true_order)
    true_pos = {g: i for i, g in enumerate(true_flat)}

    gene_to_true_chr = {}
    true_chr_boundaries = []
    current_x = 0
    cmap_true = plt.get_cmap("tab20")
    cmap_rec = plt.get_cmap("tab20b")
    true_colors = {cid: cmap_true(i % 20) for i, cid in enumerate(true_order)}
    for cid in true_order:
        start_x = current_x
        for g in true_root.get(cid, []):
            gene_to_true_chr[g] = cid
            current_x += 1
        end_x = current_x
        true_chr_boundaries.append((cid, start_x, end_x, true_colors[cid]))

    rec_genes_by_cid = rec_karyo
    rec_colors = {cid: cmap_rec(i % 20) for i, cid in enumerate(rec_order)}

    rec_flat = []
    rec_chr_boundaries = []
    current_y = 0
    for cid in rec_order:
        start_y = current_y
        for g in rec_genes_by_cid.get(cid, []):
            rec_flat.append(g)
            current_y += 1
        end_y = current_y
        rec_chr_boundaries.append((cid, start_y, end_y, rec_colors[cid]))
    rec_pos = {g: i for i, g in enumerate(rec_flat)}

    xs = []
    ys = []
    cs = []
    for g, x in true_pos.items():
        y = rec_pos.get(g)
        if y is None:
            continue
        xs.append(x)
        ys.append(y)
        cs.append(true_colors.get(gene_to_true_chr.get(g), "gray"))

    total_genes_x = len(true_flat)
    total_genes_y = len(rec_flat)

    if not xs:
        return False

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    font_family = "DejaVu Sans"
    fig, ax = plt.subplots(figsize=(15, 15))
    ax.scatter(xs, ys, c=cs, s=8, marker=".", edgecolors="none", alpha=0.6)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_title("TrueRoot", fontsize=16, pad=10, y=1.04, fontweight="bold", fontfamily=font_family)
    ax.set_ylabel("ReconstructedRoot", fontsize=16, labelpad=80, fontweight="bold", fontfamily=font_family)

    ax.set_xlim(0, total_genes_x)
    ax.set_ylim(0, total_genes_y)
    ax.invert_yaxis()
    for (cid, sx, ex, col) in true_chr_boundaries:
        ax.axvline(x=sx, color="black", linestyle="-", linewidth=0.5)
        ax.axvline(x=ex, color="black", linestyle="-", linewidth=0.5)
        bar_height = max(1, int(total_genes_y * 0.02))
        rect = mpatches.Rectangle((sx, -bar_height), ex - sx, bar_height, color=col, clip_on=False)
        ax.add_patch(rect)
        ax.text((sx + ex) / 2, -bar_height * 1.5, cid, ha="center", va="bottom", fontsize=10, fontweight="bold", fontfamily=font_family)

    for (cid, sy, ey, col) in rec_chr_boundaries:
        ax.axhline(y=sy, color="black", linestyle="-", linewidth=0.5)
        ax.axhline(y=ey, color="black", linestyle="-", linewidth=0.5)
        ax.text(-0.02, (sy + ey) / 2, cid, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=10, fontweight="bold", fontfamily=font_family, clip_on=False)

    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close(fig)
    return True

def extract_fusion_points(karyo):
    """Extract all adjacent gene pairs from a karyotype dictionary.
    
    Returns a set of normalized gene pairs (smaller gene first).
    """
    points = set()
    for genes in karyo.values():
        if not genes or len(genes) < 2:
            continue
        for g1, g2 in zip(genes, genes[1:]):
            points.add((g1, g2) if g1 <= g2 else (g2, g1))
    return points

def _normalize_pair(g1, g2):
    """Normalize gene pair order for consistent comparison."""
    return (g1, g2) if g1 <= g2 else (g2, g1)

def _find_first_index(seq, target):
    """Find first index of target in sequence, or None if not found."""
    try:
        return seq.index(target)
    except ValueError:
        return None

def _other_neighbor_in_source(source_genes, g, within_block_neighbor):
    idx = _find_first_index(source_genes, g)
    if idx is None:
        return None
    prev_g = source_genes[idx - 1] if idx - 1 >= 0 else None
    next_g = source_genes[idx + 1] if idx + 1 < len(source_genes) else None
    if within_block_neighbor is None:
        return prev_g or next_g
    if prev_g == within_block_neighbor:
        return next_g
    if next_g == within_block_neighbor:
        return prev_g
    return prev_g or next_g

def _boundary_junctions_with_triples(derived_genes, set_a, set_b):
    labels = []
    for g in derived_genes:
        if g in set_a:
            labels.append("A")
        elif g in set_b:
            labels.append("B")
        else:
            labels.append("O")

    junctions = []
    for i in range(len(derived_genes) - 1):
        la = labels[i]
        lb = labels[i + 1]
        if la not in ("A", "B") or lb not in ("A", "B") or la == lb:
            continue
        left = derived_genes[i]
        right = derived_genes[i + 1]
        left_flank = derived_genes[i - 1] if i - 1 >= 0 and labels[i - 1] == la else None
        right_flank = derived_genes[i + 2] if i + 2 < len(derived_genes) and labels[i + 2] == lb else None
        if left_flank is not None:
            triple = (left_flank, left, right)
        elif right_flank is not None:
            triple = (left, right, right_flank)
        else:
            triple = (left, left, right)
        junctions.append(
            dict(
                index=i,
                left=left,
                right=right,
                left_label=la,
                right_label=lb,
                left_flank=left_flank,
                right_flank=right_flank,
                triple=triple,
            )
        )
    return junctions

def _build_adj(genes):
    """Build {gene: set(neighbors)} adjacency dict from a gene list."""
    adj = {}
    for i in range(len(genes) - 1):
        g1, g2 = genes[i], genes[i + 1]
        adj.setdefault(g1, set()).add(g2)
        adj.setdefault(g2, set()).add(g1)
    return adj

def _analyze_rct_assuming_ancestral(ancestral_chr1, ancestral_chr2, derived_chr1, derived_chr2, a1_set=None, a2_set=None):
    """Detect RCT fusion points via adjacency scanning.
    
    Returns:
        dict with orig_pairs, trans_pairs, orig_triples, trans_triples
    """
    a1 = ancestral_chr1["genes"]
    a2 = ancestral_chr2["genes"]
    d1 = derived_chr1["genes"]
    d2 = derived_chr2["genes"]
    
    if a1_set is None: a1_set = set(a1)
    if a2_set is None: a2_set = set(a2)
    
    d1_set = set(d1)
    d2_set = set(d2)
    
    # 1. 集合差找移动基因
    moved1 = a1_set - d1_set
    moved2 = a2_set - d2_set
    
    if not moved1 or not moved2:
        return dict(trans_pairs=[], orig_pairs=[], orig_triples=[], trans_triples=[])
    
    if moved1 != (d2_set - a2_set) or moved2 != (d1_set - a1_set):
        return dict(trans_pairs=[], orig_pairs=[], orig_triples=[], trans_triples=[])
        
    trans_pairs = []
    for d in [d1, d2]:
        for i in range(len(d) - 1):
            g1, g2 = d[i], d[i+1]
            if (g1 in a1_set and g2 in a2_set) or (g1 in a2_set and g2 in a1_set):
                trans_pairs.append(_normalize_pair(g1, g2))
                
    orig_pairs = []
    for a in [a1, a2]:
        for i in range(len(a) - 1):
            g1, g2 = a[i], a[i+1]
            g1_in_d1 = g1 in d1_set
            g2_in_d1 = g2 in d1_set
            g1_in_d2 = g1 in d2_set
            g2_in_d2 = g2 in d2_set
            if (g1_in_d1 and g2_in_d2) or (g1_in_d2 and g2_in_d1):
                orig_pairs.append(_normalize_pair(g1, g2))
                
    return dict(
        trans_pairs=trans_pairs,
        orig_pairs=orig_pairs,
        orig_triples=[],
        trans_triples=[]
    )

def _outgroup_hit_details(pairs, out_adj):
    hits = []
    for p in pairs:
        if p in out_adj:
            hits.append(p)
    return hits

def normalize_id_list(*cids):
    return tuple(sorted(cids))

def root_rct_outgroup_decision(cfg, reconstructor):
    if not cfg.get("enable_root_rct_outgroup_decision", True):
        return []
    outgroup_k = getattr(reconstructor, "outgroup_karyotypes", None) or {}
    if not outgroup_k:
        reconstructor.logs.append("[Root-RCT-Outgroup-Decision] skipped: missing outgroup_karyotypes")
        return []

    out_adj = extract_fusion_points(outgroup_k)
    root = reconstructor.tree.get_tree_root()
    children = list(root.children)
    if len(children) < 2:
        reconstructor.logs.append("[Root-RCT-Outgroup-Decision] skipped: Root has <2 children")
        return []
    children.sort(key=lambda n: str(n.name))
    ln = children[0].name
    rn = children[1].name

    strict_gene_sets = bool(cfg.get("root_rct_autoscan_strict_gene_sets", True))
    require_outgroup_hit = bool(cfg.get("root_rct_autoscan_require_outgroup_hit", True))

    name_to_node = {n.name: n for n in reconstructor.tree.traverse()}

    def _ancestor_lineage_confirmed_gene_sets(node_name):
        sets_ = []
        node_t = name_to_node.get(node_name)
        if node_t is not None:
            for anc in reversed(list(node_t.iter_ancestors())):
                for genes in (reconstructor.confirmed_ancestors or {}).get(anc.name, {}).values():
                    sets_.append(frozenset(genes))
        return set(sets_)

    def _local_confirmed_lookup(node_name):
        local = (reconstructor.confirmed_ancestors or {}).get(node_name, {}) or {}
        local_sets = {}
        for cid, genes in local.items():
            local_sets[frozenset(genes)] = cid
        return local, local_sets

    def candidate_chromosomes(node_name):
        nk = reconstructor.node_karyotypes.get(node_name)
        if nk is None:
            reconstructor.logs.append(f"[RCT-Candidates] Node {node_name}: direct-only candidates=0 lineage_confirmed=0")
            return []
        ancestor_lineage_sets = _ancestor_lineage_confirmed_gene_sets(node_name)
        local_confirmed, local_confirmed_sets = _local_confirmed_lookup(node_name)
        dedup = {}
        skipped_ancestor_confirmed = 0
        skipped_ineligible = 0

        def _candidate_mode(cid, genes, attrs):
            prov = str(attrs.get("provenance", "") or "")
            telo = bool(attrs.get("telomeres", False))
            gset = frozenset(genes)
            if cid in local_confirmed or gset in local_confirmed_sets:
                return "analysis_only_confirmed"
            if "Overlap-Demoted" in prov:
                return "analysis_only_overlap"
            if telo:
                return "replaceable"
            return None

        def _mode_rank(mode):
            return {
                "analysis_only_confirmed": 3,
                "analysis_only_overlap": 2,
                "replaceable": 1,
            }.get(mode, 0)

        def _add_candidate(cid, genes, mode, source):
            gset = frozenset(genes)
            if gset in ancestor_lineage_sets and not str(mode).startswith("analysis_only"):
                return False
            prev = dedup.get(gset)
            item = (cid, list(genes), source, mode)
            if prev is None or _mode_rank(mode) > _mode_rank(prev[3]):
                dedup[gset] = item
            return True

        for cid, genes in nk.chromosomes.items():
            attrs = nk.chr_attributes.get(cid, {})
            mode = _candidate_mode(cid, genes, attrs)
            if mode is None:
                skipped_ineligible += 1
                continue
            gset = frozenset(genes)
            if gset in ancestor_lineage_sets and not str(mode).startswith("analysis_only"):
                skipped_ancestor_confirmed += 1
                continue
            _add_candidate(cid, genes, mode, node_name)

        for cid, genes in local_confirmed.items():
            _add_candidate(cid, genes, "analysis_only_confirmed", f"{node_name}[confirmed]")

        candidates = sorted(dedup.values(), key=lambda x: len(x[1]), reverse=True)
        replaceable_count = sum(1 for x in candidates if x[3] == "replaceable")
        analysis_only_count = sum(1 for x in candidates if str(x[3]).startswith("analysis_only"))
        reconstructor.logs.append(
            f"[RCT-Candidates] Node {node_name}: direct-only candidates={len(candidates)} "
            f"replaceable={replaceable_count} analysis_only={analysis_only_count} "
            f"local_confirmed={len(local_confirmed_sets)} ancestor_lineage_confirmed={len(ancestor_lineage_sets)} "
            f"skipped_ancestor_confirmed={skipped_ancestor_confirmed} skipped_ineligible={skipped_ineligible}"
        )
        return candidates

    def index_pairs(chrs, side_label):
        # We MUST keep this method, as the user stated: "RCT-Pair-Index 与这个相关吗？有必要留着不"
        # and I confirmed it's essential for the dual-side RCT detection.
        idx = {}
        chr_sets = []
        skipped_same_side_overlap = 0
        max_within_side_overlap = int(cfg.get("root_rct_pair_max_within_side_overlap", 0))
        for item in chrs:
            if len(item) >= 4:
                cid, g, node, mode = item[:4]
            elif len(item) == 3:
                cid, g, node = item
                mode = "replaceable"
            else:
                cid, g = item[:2]
                node = None
                mode = "replaceable"
            chr_sets.append((cid, g, set(g), node, mode))
        n = len(chr_sets)
        for i in range(n):
            cid1, g1, s1, node1, mode1 = chr_sets[i]
            for j in range(i + 1, n):
                cid2, g2, s2, node2, mode2 = chr_sets[j]
                overlap = s1 & s2
                if len(overlap) > max_within_side_overlap:
                    skipped_same_side_overlap += 1
                    continue
                u = s1 | s2
                if not u:
                    continue
                idx.setdefault(frozenset(u), []).append(
                    dict(chr1=cid1, genes1=g1, set1=s1, node1=node1, mode1=mode1, chr2=cid2, genes2=g2, set2=s2, node2=node2, mode2=mode2)
                )
        reconstructor.logs.append(
            f"[RCT-Pair-Index] {side_label}: unions={len(idx)} skipped_same_side_overlap={skipped_same_side_overlap} max_within_side_overlap={max_within_side_overlap}"
        )
        return idx

    results = []
    left_chrs = candidate_chromosomes(ln)
    right_chrs = candidate_chromosomes(rn)

    # --- NEW LOGIC: Shared/Nested between RCT-Candidates ---
    root_name = root.name if root.name else "Root"
    nk_root = reconstructor.node_karyotypes.get(root_name)
    left_sets = {frozenset(c[1]): c for c in left_chrs}
    right_sets = {frozenset(c[1]): c for c in right_chrs}
    
    shared_keys = set(left_sets.keys()) & set(right_sets.keys())
    for k in shared_keys:
        lc = left_sets[k]
        rc = right_sets[k]
        cid = reconstructor._pick_preferred_equivalent_id(lc[0], rc[0])
        genes = list(lc[1]) if cid == lc[0] else list(rc[1])
        if nk_root:
            nk_root.chromosomes[cid] = genes
            nk_root.chr_attributes[cid] = {"provenance": "Shared (Strict-RCT)", "telomeres": True}
        reconstructor.confirmed_ancestors.setdefault(root_name, {})[cid] = genes
        reconstructor.logs.append(f"[RCT-Shared-Strict] {lc[0]} == {rc[0]} -> {cid}")
        
    for lk, lc in left_sets.items():
        if lk in shared_keys: continue
        for rk, rc in right_sets.items():
            if rk in shared_keys: continue
            if lk.issubset(rk):
                if reconstructor._is_contiguous_block(lc[1], rc[1], strict=True):
                    cid = lc[0]
                    genes = list(lc[1])
                    if nk_root:
                        nk_root.chromosomes[cid] = genes
                        nk_root.chr_attributes[cid] = {"provenance": "Nested (Strict-RCT)", "telomeres": True}
                    reconstructor.confirmed_ancestors.setdefault(root_name, {})[cid] = genes
                    
                    residue = rk - lk
                    residue_ordered = [g for g in rc[1] if g in residue]
                    rid = reconstructor._generate_range_id(rc[0], residue_ordered)
                    idxs = [i for i, g in enumerate(rc[1]) if g in residue]
                    end_attached = False if not idxs else (min(idxs) == 0 or max(idxs) == len(rc[1]) - 1)
                    if end_attached:
                        if nk_root:
                            nk_root.chromosomes[rid] = residue_ordered
                            nk_root.chr_attributes[rid] = {"provenance": "Nested (Strict-RCT)", "telomeres": True}
                        reconstructor.confirmed_ancestors.setdefault(root_name, {})[rid] = residue_ordered
                    reconstructor.logs.append(f"[RCT-Nested-Strict] {lc[0]} in {rc[0]} -> {cid} + Residue {rid}")
            elif rk.issubset(lk):
                if reconstructor._is_contiguous_block(rc[1], lc[1], strict=True):
                    cid = rc[0]
                    genes = list(rc[1])
                    if nk_root:
                        nk_root.chromosomes[cid] = genes
                        nk_root.chr_attributes[cid] = {"provenance": "Nested (Strict-RCT)", "telomeres": True}
                    reconstructor.confirmed_ancestors.setdefault(root_name, {})[cid] = genes
                    
                    residue = lk - rk
                    residue_ordered = [g for g in lc[1] if g in residue]
                    rid = reconstructor._generate_range_id(lc[0], residue_ordered)
                    idxs = [i for i, g in enumerate(lc[1]) if g in residue]
                    end_attached = False if not idxs else (min(idxs) == 0 or max(idxs) == len(lc[1]) - 1)
                    if end_attached:
                        if nk_root:
                            nk_root.chromosomes[rid] = residue_ordered
                            nk_root.chr_attributes[rid] = {"provenance": "Nested (Strict-RCT)", "telomeres": True}
                        reconstructor.confirmed_ancestors.setdefault(root_name, {})[rid] = residue_ordered
                    reconstructor.logs.append(f"[RCT-Nested-Strict] {rc[0]} in {lc[0]} -> {cid} + Residue {rid}")

    # Re-fetch candidates after confirming
    if shared_keys or any(lk.issubset(rk) or rk.issubset(lk) for lk in left_sets for rk in right_sets):
        left_chrs = candidate_chromosomes(ln)
        right_chrs = candidate_chromosomes(rn)
    # -------------------------------------------------------

    if len(left_chrs) < 2 or len(right_chrs) < 2:
        reconstructor.logs.append("[Root-RCT-Outgroup-Decision] skipped: insufficient direct telomeric non-confirmed chromosomes in Root children")
        return []

    L = index_pairs(left_chrs, ln)
    R = index_pairs(right_chrs, rn)
    common_unions = set(L.keys()) & set(R.keys())
    if not common_unions:
        reconstructor.logs.append("[Root-RCT-Outgroup-Decision] no matching union gene-sets across direct Root children")
        return []

    left_union_sizes = sorted({len(u) for u in L.keys()}, reverse=True)
    right_union_sizes = sorted({len(u) for u in R.keys()}, reverse=True)

    occupied_genes = set()
    for genes in reconstructor.confirmed_ancestors.get(root_name, {}).values():
        occupied_genes.update(genes)

    seen_events = set()
    seen_ancestral_sides = set()
    filter_counts = {"strict_or_intersection": 0, "require_outgroup_hit": 0, "tie": 0, "ancestral_dedup": 0, "event_dedup": 0, "emitted": 0}
    for u in sorted(common_unions, key=lambda s: len(s), reverse=True):
        if u.issubset(occupied_genes):
            continue
        left_pairs = L[u]
        right_pairs = R[u]
        for lp in left_pairs:
            l1 = dict(chr=lp["chr1"], genes=lp["genes1"])
            l2 = dict(chr=lp["chr2"], genes=lp["genes2"])
            lp_s1 = lp["set1"]
            lp_s2 = lp["set2"]
            for rp in right_pairs:
                r1 = dict(chr=rp["chr1"], genes=rp["genes1"])
                r2 = dict(chr=rp["chr2"], genes=rp["genes2"])
                rp_s1 = rp["set1"]
                rp_s2 = rp["set2"]
                if strict_gene_sets:
                    if not (lp_s1 & rp_s1 and lp_s1 & rp_s2):
                        continue
                    if not (lp_s2 & rp_s1 and lp_s2 & rp_s2):
                        continue
                else:
                    if not ((lp_s1 & rp_s1 and lp_s2 & rp_s2) or (lp_s1 & rp_s2 and lp_s2 & rp_s1)):
                        filter_counts["strict_or_intersection"] += 1
                        continue

                a1 = _analyze_rct_assuming_ancestral(l1, l2, r1, r2, a1_set=lp_s1, a2_set=lp_s2)
                b1 = _analyze_rct_assuming_ancestral(r1, r2, l1, l2, a1_set=rp_s1, a2_set=rp_s2)
                
                if not (b1["orig_pairs"] or b1["trans_pairs"]):
                    b1 = dict(
                        orig_pairs=list(a1["trans_pairs"]),
                        trans_pairs=list(a1["orig_pairs"]),
                        orig_triples=list(a1["trans_triples"]),
                        trans_triples=list(a1["orig_triples"]),
                    )
                if not (a1["orig_pairs"] or a1["trans_pairs"]):
                    a1 = dict(
                        orig_pairs=list(b1["trans_pairs"]),
                        trans_pairs=list(b1["orig_pairs"]),
                        orig_triples=list(b1["trans_triples"]),
                        trans_triples=list(b1["orig_triples"]),
                    )
                a_orig_hits = _outgroup_hit_details(a1["orig_pairs"], out_adj)
                b_orig_hits = _outgroup_hit_details(b1["orig_pairs"], out_adj)

                if require_outgroup_hit and not (a_orig_hits or b_orig_hits):
                    reconstructor.logs.append(f"[Root-RCT-Outgroup-Decision] Filtered by outgroup: {l1['chr']} + {l2['chr']} vs {r1['chr']} + {r2['chr']}")
                    reconstructor.logs.append(f"  a1 orig_pairs: {a1['orig_pairs']}")
                    reconstructor.logs.append(f"  b1 orig_pairs: {b1['orig_pairs']}")
                    reconstructor.logs.append(f"  a_orig_hits: {a_orig_hits}")
                    reconstructor.logs.append(f"  b_orig_hits: {b_orig_hits}")
                    filter_counts["require_outgroup_hit"] += 1
                    continue

                a_score = len(a_orig_hits)
                b_score = len(b_orig_hits)

                if a_score > b_score:
                    verdict = "left_is_ancestral"
                elif b_score > a_score:
                    verdict = "right_is_ancestral"
                else:
                    reconstructor.logs.append(
                        f"[Root-RCT-Outgroup-Decision] tie skipped: {ln}({l1['chr']}, {l2['chr']}) vs {rn}({r1['chr']}, {r2['chr']})"
                    )
                    filter_counts["tie"] += 1
                    continue

                cid_l1, cid_l2 = normalize_id_list(l1["chr"], l2["chr"])
                cid_r1, cid_r2 = normalize_id_list(r1["chr"], r2["chr"])
                left_analysis_only = str(lp.get("mode1", "")).startswith("analysis_only") or str(lp.get("mode2", "")).startswith("analysis_only")
                right_analysis_only = str(rp.get("mode1", "")).startswith("analysis_only") or str(rp.get("mode2", "")).startswith("analysis_only")
                winner_analysis_only = left_analysis_only if verdict == "left_is_ancestral" else right_analysis_only
                ancestral_side_key = f"{ln}({cid_l1},{cid_l2})" if verdict == "left_is_ancestral" else f"{rn}({cid_r1},{cid_r2})"
                if ancestral_side_key in seen_ancestral_sides:
                    filter_counts["ancestral_dedup"] += 1
                    continue
                seen_ancestral_sides.add(ancestral_side_key)

                cid = f"{ln}({cid_l1},{cid_l2})__{rn}({cid_r1},{cid_r2})"
                if cid in seen_events:
                    filter_counts["event_dedup"] += 1
                    continue
                seen_events.add(cid)

                winner_chroms = (
                    f"{l1['chr']}, {l2['chr']}" if verdict == "left_is_ancestral" else f"{r1['chr']}, {r2['chr']}"
                )
                reconstructor.logs.append(f"[Root-RCT-Outgroup-Decision] {cid}")
                reconstructor.logs.append(f"  RCT Event: {ln}({l1['chr']}, {l2['chr']}) vs {rn}({r1['chr']}, {r2['chr']})")
                reconstructor.logs.append(f"  H1 proto={ln}({l1['chr']}, {l2['chr']}) orig_pairs={sorted(a1['orig_pairs'])} trans_pairs={sorted(a1['trans_pairs'])} outgroup_hit={sorted(a_orig_hits)} match={len(a_orig_hits)}/{len(a1['orig_pairs'])} mode={'analysis-only' if left_analysis_only else 'replaceable'}")
                reconstructor.logs.append(f"  H2 proto={rn}({r1['chr']}, {r2['chr']}) orig_pairs={sorted(b1['orig_pairs'])} trans_pairs={sorted(b1['trans_pairs'])} outgroup_hit={sorted(b_orig_hits)} match={len(b_orig_hits)}/{len(b1['orig_pairs'])} mode={'analysis-only' if right_analysis_only else 'replaceable'}")
                reconstructor.logs.append(f"  Outgroup(adj)={sorted(list(set(a_orig_hits + b_orig_hits)))}")
                if winner_analysis_only:
                    reconstructor.logs.append(f"  Conclusion: {winner_chroms} are Ancestor [analysis-only, no replacement]")
                else:
                    reconstructor.logs.append(f"  Conclusion: {winner_chroms} are Ancestor")
                reconstructor.logs.append("")

                filter_counts["emitted"] += 1
                results.append(
                    dict(
                        ComparisonID=cid,
                        LeftNode=ln,
                        LeftChr1=l1["chr"],
                        LeftChr2=l2["chr"],
                        RightNode=rn,
                        RightChr1=r1["chr"],
                        RightChr2=r2["chr"],
                        LeftScore=a_score,
                        RightScore=b_score,
                        Verdict=verdict,
                        LeftOrigPairs=";".join([f"{x}|{y}" for x, y in sorted(a1["orig_pairs"])]),
                        LeftTransPairs=";".join([f"{x}|{y}" for x, y in sorted(a1["trans_pairs"])]),
                        LeftOutgroupHit=";".join([f"{x}|{y}" for x, y in sorted(a_orig_hits)]),
                        RightOrigPairs=";".join([f"{x}|{y}" for x, y in sorted(b1["orig_pairs"])]),
                        RightTransPairs=";".join([f"{x}|{y}" for x, y in sorted(b1["trans_pairs"])]),
                        RightOutgroupHit=";".join([f"{x}|{y}" for x, y in sorted(b_orig_hits)]),
                        LeftAnalysisOnly=left_analysis_only,
                        RightAnalysisOnly=right_analysis_only,
                        WinnerAnalysisOnly=winner_analysis_only,
                    )
                )
    return results

# -----------------------------
# 3. RECONSTRUCTION ENGINE
# -----------------------------
class NodeKaryotype:
    def __init__(self, name, chromosomes=None, is_leaf=False):
        self.name = name
        self.chromosomes = chromosomes if chromosomes else {}
        self.chr_attributes = {}
        self.is_leaf = is_leaf
        self.inference_log = []
        
        for cid in self.chromosomes:
            self.chr_attributes[cid] = {"provenance": "Unknown", "telomeres": False}

    def add_log(self, msg):
        self.inference_log.append(msg)

class AncestorReconstructor:
    def __init__(self, tree, species_karyotypes, gene_position_map, cfg=None):
        self.cfg = cfg if cfg is not None else CONFIG
        self.tree = tree
        self.species_karyotypes = species_karyotypes
        self.gene_position_map = gene_position_map 
        self.node_karyotypes = {} 
        self.final_ancestors = [] 
        self.gene_to_ancestor = {}
        self.logs = []
        self.confirmed_ancestors = {}
        self.resolved_chromosomes = {}
        self.residual_chromosomes = {}
        self.pending_ancestors = {}
        self.pending_gene_sets_by_node = {}
        self._residual_pool_build_count = 0
        self.frozen_confirmed_ancestors = {}
        self.frozen_confirmed_gene_sets = {}
        self._short_circuit_after_strict = False
        
    def _purge_residual_nested_supersets(self):
        for node_name, nk in self.node_karyotypes.items():
            to_remove = [cid for cid, attrs in nk.chr_attributes.items() if attrs.get("provenance") == "Residual (Nested-Superset)"]
            for cid in to_remove:
                nk.chromosomes.pop(cid, None)
                nk.chr_attributes.pop(cid, None)
                if node_name in self.confirmed_ancestors:
                    self.confirmed_ancestors[node_name].pop(cid, None)
                    if not self.confirmed_ancestors[node_name]:
                        self.confirmed_ancestors.pop(node_name, None)

    def _get_node_by_name(self, node_name):
        for node in self.tree.traverse():
            if node.name == node_name:
                return node
        return None

    def _get_applicable_confirmed_gene_sets(self, node_name, include_self=True):
        node = self._get_node_by_name(node_name)
        if node is None:
            return set()
        applicable = set()
        if include_self:
            for genes in self.confirmed_ancestors.get(node_name, {}).values():
                applicable.add(frozenset(genes))
        for ancestor in node.iter_ancestors():
            for genes in self.confirmed_ancestors.get(ancestor.name, {}).values():
                applicable.add(frozenset(genes))
        return applicable

    def _rebuild_residual_pools(self, reason=""):
        self._init_residual_pools()
        self._residual_pool_build_count += 1

    def _snapshot_frozen_confirmed_ancestors(self, reason=""):
        frozen = {}
        frozen_sets = {}
        total = 0
        for node_name, node_map in (self.confirmed_ancestors or {}).items():
            if not node_map:
                continue
            frozen[node_name] = {cid: list(genes) for cid, genes in node_map.items()}
            frozen_sets[node_name] = {frozenset(genes): cid for cid, genes in node_map.items()}
            total += len(node_map)
        self.frozen_confirmed_ancestors = frozen
        self.frozen_confirmed_gene_sets = frozen_sets
        if reason:
            self.logs.append(f"[Freeze-Confirmed] reason={reason} nodes={len(frozen)} total={total}")

    def _get_frozen_cid_for_genes(self, node_name, genes):
        return (self.frozen_confirmed_gene_sets.get(node_name, {}) or {}).get(frozenset(genes))

    def _is_frozen_confirmed(self, node_name, cid=None, genes=None):
        node_map = self.frozen_confirmed_ancestors.get(node_name, {}) or {}
        if cid is not None and cid in node_map:
            if genes is None:
                return True
            return set(node_map[cid]) == set(genes)
        if genes is not None:
            return frozenset(genes) in (self.frozen_confirmed_gene_sets.get(node_name, {}) or {})
        return False

    def _root_gene_repertoire_complete(self):
        root = self.tree.get_tree_root()
        if root is None:
            self.logs.append("[Root-Gene-Repertoire-Check] skipped: tree has no root")
            return False

        root_name = root.name
        root_confirmed = self.confirmed_ancestors.get(root_name, {})
        root_genes = set()
        for genes in root_confirmed.values():
            root_genes.update(genes)

        modern_genes = set()
        if self.species_karyotypes:
            sp = list(self.species_karyotypes.keys())[0]
            for genes in self.species_karyotypes[sp].values():
                modern_genes.update(genes)

        root_len = len(root_genes)
        modern_len = len(modern_genes)
        if root_len != modern_len:
            self.logs.append(
                f"[Root-Gene-Repertoire-Check] length_mismatch root={root_len} modern={modern_len}"
            )
            return False

        if root_genes != modern_genes:
            missing = sorted(modern_genes - root_genes)
            extra = sorted(root_genes - modern_genes)
            missing_preview = ",".join(map(str, missing[:10])) if missing else ""
            extra_preview = ",".join(map(str, extra[:10])) if extra else ""
            self.logs.append(
                f"[Root-Gene-Repertoire-Check] set_mismatch missing={len(missing)} extra={len(extra)} "
                f"missing_preview={missing_preview} extra_preview={extra_preview}"
            )
            return False

        self.logs.append(
            f"[Root-Gene-Repertoire-Check] complete: root={root_len} modern={modern_len}"
        )
        return True

    def _finalize_root_pending_ancestors(self):
        if not self.cfg.get("enable_root_final_pending_promotion", True):
            return
        if not self.pending_ancestors:
            return
        min_count = int(self.cfg.get("root_final_pending_min_count", 2))
        root = self.tree.get_tree_root()
        root_name = root.name
        root_label = "Root" if root_name else "Root"
        nk_root = self.node_karyotypes.get(root_name)
        if nk_root is None:
            return

        root_children = [child.name for child in root.children if not child.is_leaf()]
        if not root_children:
            return

        derived_nodes = set()
        if self.cfg.get("enable_root_rct_outgroup_decision", True):
            for row in getattr(self, "root_rct_outgroup_decisions", []) or []:
                left = row.get("LeftNode")
                right = row.get("RightNode")
                verdict = row.get("Verdict")
                if verdict == "left_is_ancestral" and right:
                    derived_nodes.add(right)
                if verdict == "right_is_ancestral" and left:
                    derived_nodes.add(left)

        existing_sets = {frozenset(genes) for genes in self.confirmed_ancestors.get(root_name, {}).values()}
        occupied = set()
        for genes in self.confirmed_ancestors.get(root_name, {}).values():
            occupied.update(genes)

        root_ids_by_set = {}
        for cid, genes in nk_root.chromosomes.items():
            root_ids_by_set.setdefault(frozenset(genes), []).append(cid)

        promoted = 0
        for sp_prefix, d in self.pending_ancestors.items():
            for key, entry in d.items():
                if key in existing_sets:
                    continue
                if entry.get("count", 0) < min_count:
                    continue
                if derived_nodes:
                    ev = entry.get("evidence", [])
                    if ev and all(any(f"{dn}:" in e for dn in derived_nodes) for e in ev):
                        ev_note = "; ".join(ev[:2])
                        self.logs.append(f"[Final-Pending-Skip-Derived] {root_label} sp={sp_prefix} count={entry.get('count',0)} evidence={ev_note}")
                        continue
                
                genes_set = set(key)
                if not genes_set:
                    continue
                if genes_set & occupied:
                    continue
                
                validation_by_child = {}
                for child_name in root_children:
                    child_confirmed = self.confirmed_ancestors.get(child_name, {})
                    child_has_match = False
                    for cid, cgenes in child_confirmed.items():
                        if set(cgenes) == genes_set:
                            child_has_match = True
                            validation_by_child[child_name] = f"confirmed:{cid}"
                            break
                    
                    if not child_has_match:
                        child_residuals = self.residual_chromosomes.get(child_name, {})
                        for cid, cgenes in child_residuals.items():
                            if set(cgenes) == genes_set:
                                child_has_match = True
                                validation_by_child[child_name] = f"residual:{cid}"
                                break
                            if genes_set.issubset(set(cgenes)):
                                child_has_match = True
                                validation_by_child[child_name] = f"residual_subset:{cid}"
                                break
                    
                    if not child_has_match:
                        validation_by_child[child_name] = None
                
                unvalidated_children = [name for name, val in validation_by_child.items() if val is None]
                if unvalidated_children:
                    ev_note = "; ".join(entry.get("evidence", [])[:2]) if entry.get("evidence") else "N/A"
                    self.logs.append(f"[Final-Pending-Skip-Unvalidated] {root_label} sp={sp_prefix} genes={len(genes_set)} unvalidated_by={', '.join(unvalidated_children)} evidence={ev_note}")
                    continue
                existing_root_ids = root_ids_by_set.get(key, [])
                promoted_id = None
                for cid in existing_root_ids:
                    attrs = nk_root.chr_attributes.get(cid, {})
                    if not attrs.get("telomeres", False):
                        promoted_id = cid
                        break

                if promoted_id is not None:
                    nk_root.chr_attributes[promoted_id] = {"provenance": "Reconciled (Strict-RCT)", "telomeres": True}
                    self.confirmed_ancestors.setdefault(root_name, {})[promoted_id] = list(nk_root.chromosomes.get(promoted_id, []))
                else:
                    genes = list(entry.get("genes") or genes_set)
                    if sp_prefix in self.gene_position_map:
                        sp_map = self.gene_position_map[sp_prefix]
                        genes = sorted(genes_set, key=lambda g: sp_map.get(g, ("", 10**12))[1])
                    new_id = self._generate_range_id(f"{sp_prefix}(pending_final)", genes)
                    if new_id in nk_root.chromosomes and frozenset(nk_root.chromosomes.get(new_id, [])) != key:
                        new_id = f"{new_id}_final"
                    nk_root.chromosomes[new_id] = list(genes)
                    nk_root.chr_attributes[new_id] = {"provenance": "Reconciled (Strict-RCT)", "telomeres": True}
                    self.confirmed_ancestors.setdefault(root_name, {})[new_id] = list(genes)
                    root_ids_by_set.setdefault(key, []).append(new_id)
                    promoted_id = new_id
                existing_sets.add(key)
                occupied.update(genes_set)
                ev = entry.get("evidence", [])
                ev_note = "; ".join(ev[:2])
                self.logs.append(f"[Final-Pending-Promote] {root_label}:{promoted_id} sp={sp_prefix} count={entry.get('count',0)} evidence={ev_note}")
                promoted += 1

        if promoted:
            self.logs.append(f"[Final-Pending-Promote-Summary] {root_label} promoted={promoted}")

    def _promote_unmatched_ancestors_via_outgroup(self):
        """
        Identify confirmed ancestral chromosomes in Anc1 (or Anc2) that split exactly into 2 contiguous blocks 
        in the sibling node's confirmed chromosomes, OR map to fragments in extant sibling species.
        Check the breakpoint connection against the outgroup:
        - If outgroup supports the connection -> The single chromosome is the true Ancestor (RCT occurred).
        - If outgroup does NOT support it -> The two separate blocks are the true Ancestors (EEJ occurred).
        """
        root = self.tree.get_tree_root()
        children = [child for child in root.children if not child.is_leaf()]
        if len(children) < 2:
            return
            
        outgroup_k = getattr(self, "outgroup_karyotypes", None) or {}
        if not outgroup_k:
            return
        out_adj = extract_fusion_points(outgroup_k)

        root_name = root.name
        nk_root = self.node_karyotypes.get(root_name)
        if nk_root is None:
            return

        def _get_extant_descendants(node):
            return [n for n in node.traverse() if n.is_leaf()]

        promoted_count = 0
        
        for main_node, sibling_node in [(children[0], children[1]), (children[1], children[0])]:
            main_confirmed = dict(self.confirmed_ancestors.get(main_node.name, {}))
            # Also include telomeres=True residue chromosomes from the node karyotype:
            # these are chromosomes that exist in the node's inferred karyotype but were
            # not confirmed via Shared/Nested (e.g., Sp14(3):375-891 at Anc2 is a Residue
            # with telomeres=True that Logic 1 should be able to promote).
            main_nk = self.node_karyotypes.get(main_node.name)
            if main_nk:
                confirmed_gsets = {frozenset(g) for g in main_confirmed.values()}
                for cid, genes in main_nk.chromosomes.items():
                    if cid not in main_confirmed:
                        attr = main_nk.chr_attributes.get(cid, {})
                        if attr.get("telomeres", False):
                            gset_cid = frozenset(genes)
                            if gset_cid not in confirmed_gsets:
                                confirmed_gsets.add(gset_cid)
                                main_confirmed[cid] = genes

            sibling_confirmed = self.confirmed_ancestors.get(sibling_node.name, {})

            root_sets = {frozenset(genes) for genes in self.confirmed_ancestors.get(root_name, {}).values()}

            extant_siblings = _get_extant_descendants(sibling_node)

            # ---------------------------------------------------------
            # LOGIC 1: Single Chromosome Validation
            # ---------------------------------------------------------
            for cid, genes in main_confirmed.items():
                gset = frozenset(genes)
                if any(gset & rset for rset in root_sets):
                    continue

                # Skip if the sibling already has this chromosome (would be handled by Shared-Strict)
                shares_genes_with_sibling = False
                for scid, s_genes in sibling_confirmed.items():
                    if gset & frozenset(s_genes):
                        shares_genes_with_sibling = True
                        break
                if shares_genes_with_sibling:
                    continue

                # Reject multi-family chromosomes
                def _gene_family_prefix(g):
                    s = str(g)
                    i = 0
                    while i < len(s) and s[i].isalpha():
                        i += 1
                    return s[:i] if i > 0 else ""

                _fams = {_gene_family_prefix(g) for g in genes}
                if len(_fams) > 1:
                    continue

                valid_breakpoint = None
                outgroup_rejected_breakpoint = None

                # Check how this chromosome distributes in sibling's extant descendants
                for extant in extant_siblings:
                    extant_nk = self.node_karyotypes.get(extant.name)
                    if not extant_nk:
                        continue

                    overlapping_chrs = []
                    for ecid, egenes in extant_nk.chromosomes.items():
                        eset = frozenset(egenes)
                        intersection = eset & gset
                        if intersection:
                            overlapping_chrs.append((ecid, egenes, intersection))

                    # User requirement: 断裂为两部分的要求是：这两部分基因数量总和与染色体的基因数量相等
                    if len(overlapping_chrs) == 2:
                        p1_genes = overlapping_chrs[0][2]
                        p2_genes = overlapping_chrs[1][2]

                        if len(p1_genes) + len(p2_genes) == len(genes) and (p1_genes | p2_genes) == gset:
                            cross_points = []
                            gene_to_piece = {}
                            for idx, (_, _, inter) in enumerate(overlapping_chrs):
                                for g in inter:
                                    gene_to_piece[g] = idx
                                    
                            for i in range(len(genes) - 1):
                                g1, g2 = genes[i], genes[i+1]
                                if g1 in gene_to_piece and g2 in gene_to_piece:
                                    if gene_to_piece[g1] != gene_to_piece[g2]:
                                        cross_points.append(_normalize_pair(g1, g2))
                                        
                            if len(cross_points) == 1:
                                bp = cross_points[0]
                                if bp in out_adj:
                                    valid_breakpoint = bp
                                    break
                                else:
                                    outgroup_rejected_breakpoint = bp
                            
                if valid_breakpoint:
                    self.logs.append(f"[Isolated-Outgroup-Validated] {cid} from {main_node.name} is validated by outgroup breakpoint {valid_breakpoint} mapping to fragments in sibling lineage. Promoted to Root.")
                    new_id = self._generate_range_id(cid, genes)
                    if nk_root:
                        nk_root.chromosomes[new_id] = list(genes)
                        nk_root.chr_attributes[new_id] = {"provenance": "Isolated-Outgroup-Validated", "telomeres": True}
                    self.confirmed_ancestors.setdefault(root_name, {})[new_id] = list(genes)
                    root_sets.add(gset)
                    promoted_count += 1

            # ---------------------------------------------------------
            # LOGIC 2 & 3: Cross-Chromosome Fusion Validation & RCT Fallback
            # ---------------------------------------------------------
            for cid, genes in main_confirmed.items():
                gset = frozenset(genes)
                if gset in root_sets:
                    continue

                overlapping_chrs = []
                for scid, sgenes in sibling_confirmed.items():
                    eset = frozenset(sgenes)
                    intersection = eset & gset
                    if intersection:
                        overlapping_chrs.append((scid, sgenes, intersection))
                        
                if len(overlapping_chrs) == 2:
                    p1_genes = overlapping_chrs[0][2]
                    p2_genes = overlapping_chrs[1][2]
                    
                    if len(p1_genes) + len(p2_genes) == len(genes) and (p1_genes | p2_genes) == gset:
                        cross_points = []
                        gene_to_piece = {}
                        for idx, (_, _, inter) in enumerate(overlapping_chrs):
                            for g in inter:
                                gene_to_piece[g] = idx
                                
                        for i in range(len(genes) - 1):
                            g1, g2 = genes[i], genes[i+1]
                            if g1 in gene_to_piece and g2 in gene_to_piece:
                                if gene_to_piece[g1] != gene_to_piece[g2]:
                                    cross_points.append(_normalize_pair(g1, g2))
                                    
                        if len(cross_points) == 1:
                            bp = cross_points[0]
                            if bp in out_adj:
                                if not any(gset & rset for rset in root_sets):
                                    self.logs.append(f"[Fusion-Outgroup-Validated] {cid} from {main_node.name} splits exactly into 2 confirmed ancestors in {sibling_node.name}. Breakpoint {bp} supported by outgroup. Promoted to Root.")
                                    new_id = self._generate_range_id(cid, genes)
                                    if nk_root:
                                        nk_root.chromosomes[new_id] = list(genes)
                                        nk_root.chr_attributes[new_id] = {"provenance": "Fusion-Outgroup-Validated", "telomeres": True}
                                    self.confirmed_ancestors.setdefault(root_name, {})[new_id] = list(genes)
                                    root_sets.add(gset)
                                    promoted_count += 1
                            else:
                                # Logic 3: RCT Fallback for Fusion
                                found_rct_partner = False
                                for other_cid, other_genes in main_confirmed.items():
                                    if other_cid == cid:
                                        continue
                                    
                                    u = gset | frozenset(other_genes)
                                    current_connections = set()
                                    for g_list in [genes, other_genes]:
                                        for k in range(len(g_list) - 1):
                                            current_connections.add(_normalize_pair(g_list[k], g_list[k+1]))
                                            
                                    if bp not in current_connections:
                                        continue
                                        
                                    supported_connections = set()
                                    for g1 in gset:
                                        for g2 in other_genes:
                                            pair = _normalize_pair(g1, g2)
                                            if pair not in current_connections and pair in out_adj:
                                                supported_connections.add(pair)
                                                
                                    if len(supported_connections) == 2:
                                        # Outgroup supports the CROSS connections, which exist in the sibling node.
                                        # This means the current node's chromosomes (cid and other_cid) are the DERIVED state (fused).
                                        # The true ancestors are the blocks in the sibling node that form the cross connections.
                                        # Find which sibling confirmed chromosomes contain these genes
                                        
                                        sibling_ancestors_to_promote = []
                                        for scid, sgenes in sibling_confirmed.items():
                                            sset = frozenset(sgenes)
                                            # If this sibling chromosome is part of the union of the 2 derived chromosomes
                                            if sset.issubset(u):
                                                # Check if it contains any of the supported cross-connections
                                                contains_supported = False
                                                for p in supported_connections:
                                                    if p[0] in sset and p[1] in sset:
                                                        for k in range(len(sgenes) - 1):
                                                            if _normalize_pair(sgenes[k], sgenes[k+1]) == p:
                                                                contains_supported = True
                                                                break
                                                if contains_supported:
                                                    sibling_ancestors_to_promote.append((scid, sgenes))
                                                    
                                        if sibling_ancestors_to_promote:
                                            promoted_ids = [c for c, _ in sibling_ancestors_to_promote]
                                            self.logs.append(f"[Fusion-Outgroup-Rejected-RCT] {cid} breakpoint {bp} rejected, but forms RCT with {other_cid}. True ancestors are in {sibling_node.name}: {promoted_ids}.")
                                            for promote_cid, promote_genes in sibling_ancestors_to_promote:
                                                p_gset = frozenset(promote_genes)
                                                if not any(p_gset & rset for rset in root_sets):
                                                    new_id = self._generate_range_id(promote_cid, promote_genes)
                                                    if nk_root:
                                                        nk_root.chromosomes[new_id] = list(promote_genes)
                                                        nk_root.chr_attributes[new_id] = {"provenance": "Fusion-Outgroup-Rejected-RCT", "telomeres": True}
                                                    self.confirmed_ancestors.setdefault(root_name, {})[new_id] = list(promote_genes)
                                                    root_sets.add(p_gset)
                                                    promoted_count += 1
                                            found_rct_partner = True
                                            break
                                        
                                if found_rct_partner:
                                    continue

            # ---------------------------------------------------------
            # LOGIC 2 & 3: Cross-Chromosome Fusion Validation & RCT Fallback
            # ---------------------------------------------------------
            for cid, genes in main_confirmed.items():
                gset = frozenset(genes)
                if gset in root_sets:
                    continue

                overlapping_chrs = []
                for scid, sgenes in sibling_confirmed.items():
                    eset = frozenset(sgenes)
                    intersection = eset & gset
                    if intersection:
                        overlapping_chrs.append((scid, sgenes, intersection))
                        
                if len(overlapping_chrs) == 2:
                    p1_genes = overlapping_chrs[0][2]
                    p2_genes = overlapping_chrs[1][2]
                    
                    if len(p1_genes) + len(p2_genes) == len(genes) and (p1_genes | p2_genes) == gset:
                        cross_points = []
                        gene_to_piece = {}
                        for idx, (_, _, inter) in enumerate(overlapping_chrs):
                            for g in inter:
                                gene_to_piece[g] = idx
                                
                        for i in range(len(genes) - 1):
                            g1, g2 = genes[i], genes[i+1]
                            if g1 in gene_to_piece and g2 in gene_to_piece:
                                if gene_to_piece[g1] != gene_to_piece[g2]:
                                    cross_points.append(_normalize_pair(g1, g2))
                                    
                        if len(cross_points) == 1:
                            bp = cross_points[0]
                            if bp in out_adj:
                                if not any(gset & rset for rset in root_sets):
                                    self.logs.append(f"[Fusion-Outgroup-Validated] {cid} from {main_node.name} splits exactly into 2 confirmed ancestors in {sibling_node.name}. Breakpoint {bp} supported by outgroup. Promoted to Root.")
                                    new_id = self._generate_range_id(cid, genes)
                                    if nk_root:
                                        nk_root.chromosomes[new_id] = list(genes)
                                        nk_root.chr_attributes[new_id] = {"provenance": "Fusion-Outgroup-Validated", "telomeres": True}
                                    self.confirmed_ancestors.setdefault(root_name, {})[new_id] = list(genes)
                                    root_sets.add(gset)
                                    promoted_count += 1
                            else:
                                # Logic 3: RCT Fallback for Fusion
                                found_rct_partner = False
                                for other_cid, other_genes in main_confirmed.items():
                                    if other_cid == cid:
                                        continue
                                    
                                    current_connections = set()
                                    for g_list in [genes, other_genes]:
                                        for k in range(len(g_list) - 1):
                                            current_connections.add(_normalize_pair(g_list[k], g_list[k+1]))
                                            
                                    if bp not in current_connections:
                                        continue
                                        
                                    supported_connections = set()
                                    for g1 in gset:
                                        for g2 in other_genes:
                                            pair = _normalize_pair(g1, g2)
                                            if pair not in current_connections and pair in out_adj:
                                                supported_connections.add(pair)
                                                
                                    if len(supported_connections) == 2:
                                        sibling_ancestors_to_promote = []
                                        for scid_sib, sgenes_sib in sibling_confirmed.items():
                                            sset_sib = frozenset(sgenes_sib)
                                            if sset_sib.issubset(gset | frozenset(other_genes)):
                                                contains_supported = False
                                                for p in supported_connections:
                                                    if p[0] in sset_sib and p[1] in sset_sib:
                                                        for k in range(len(sgenes_sib) - 1):
                                                            if _normalize_pair(sgenes_sib[k], sgenes_sib[k+1]) == p:
                                                                contains_supported = True
                                                                break
                                                if contains_supported:
                                                    sibling_ancestors_to_promote.append((scid_sib, sgenes_sib))
                                                    
                                        if sibling_ancestors_to_promote:
                                            promoted_ids = [c for c, _ in sibling_ancestors_to_promote]
                                            self.logs.append(f"[Fusion-Outgroup-Rejected-RCT] {cid} breakpoint {bp} rejected, but forms RCT with {other_cid}. True ancestors are in {sibling_node.name}: {promoted_ids}.")
                                            for promote_cid, promote_genes in sibling_ancestors_to_promote:
                                                p_gset = frozenset(promote_genes)
                                                if not any(p_gset & rset for rset in root_sets):
                                                    new_id = self._generate_range_id(promote_cid, promote_genes)
                                                    if nk_root:
                                                        nk_root.chromosomes[new_id] = list(promote_genes)
                                                        nk_root.chr_attributes[new_id] = {"provenance": "Fusion-Outgroup-Rejected-RCT", "telomeres": True}
                                                    self.confirmed_ancestors.setdefault(root_name, {})[new_id] = list(promote_genes)
                                                    root_sets.add(p_gset)
                                                    promoted_count += 1
                                            found_rct_partner = True
                                            break
                                        
                                if found_rct_partner:
                                    continue

            # ---------------------------------------------------------
            # LOGIC 4: Outgroup-Arbitrated Nested Validation
            # ---------------------------------------------------------
            for cid_sub, genes_sub in main_confirmed.items():
                set_sub = frozenset(genes_sub)
                if any(set_sub & rset for rset in root_sets):
                    continue

                for cid_super, genes_super in sibling_confirmed.items():
                    set_super = frozenset(genes_super)
                    if set_super in root_sets:
                        continue
                        
                    # Check if it's a strict subset
                    if set_sub < set_super:
                        blocks = []
                        current_block = []
                        for g in genes_super:
                            if g in set_sub:
                                current_block.append(g)
                            else:
                                if current_block:
                                    blocks.append(current_block)
                                    current_block = []
                        if current_block:
                            blocks.append(current_block)
                            
                        # If it's split into exactly 2 blocks
                        if len(blocks) == 2:
                            b1_set = frozenset(blocks[0])
                            b2_set = frozenset(blocks[1])
                            
                            cross_points = []
                            for i in range(len(genes_sub) - 1):
                                g1, g2 = genes_sub[i], genes_sub[i+1]
                                if (g1 in b1_set and g2 in b2_set) or (g1 in b2_set and g2 in b1_set):
                                    cross_points.append(_normalize_pair(g1, g2))
                                    
                            if len(cross_points) == 1:
                                bp = cross_points[0]
                                if bp in out_adj:
                                    self.logs.append(f"[Outgroup-Arbitrated-Nested] {cid_sub} ({main_node.name}) is a true subset of {cid_super} ({sibling_node.name}), split in 2 blocks. Outgroup supports {bp}. {cid_sub} is the true Ancestor!")
                                    
                                    new_id = self._generate_range_id(cid_sub, genes_sub)
                                    if nk_root:
                                        nk_root.chromosomes[new_id] = list(genes_sub)
                                        nk_root.chr_attributes[new_id] = {"provenance": "Outgroup-Arbitrated-Nested", "telomeres": True}
                                    self.confirmed_ancestors.setdefault(root_name, {})[new_id] = list(genes_sub)
                                    root_sets.add(set_sub)
                                    promoted_count += 1
                                    
                                    residue_genes = [g for g in genes_super if g not in set_sub]
                                    if residue_genes:
                                        res_id = self._generate_range_id(f"{cid_super}_residue", residue_genes)
                                        if sibling_node.name not in self.residual_chromosomes:
                                            self.residual_chromosomes[sibling_node.name] = {}
                                        self.residual_chromosomes[sibling_node.name][res_id] = residue_genes
                                        self.logs.append(f"  -> Generated residue {res_id} in {sibling_node.name} from the remaining parts of {cid_super}.")

            # ---------------------------------------------------------
            # LOGIC 5: RCT Residual Reconstruction via Outgroup Adjacency
            # ---------------------------------------------------------
            # For each confirmed Ci at main_node:
            #   1. Compute Ci_residual = genes in Ci not yet root-confirmed
            #   2. For each Cj in the same node's full karyotype:
            #      a. Compute Cj_residual similarly (disjoint from Ci_residual)
            #      b. Find out_adj pairs crossing Ci_residual ↔ Cj_residual
            #      c. Exactly 1 such cross-connection → merge into a new Root ancestor
            # ---------------------------------------------------------
            main_nk = self.node_karyotypes.get(main_node.name)
            if main_nk is not None:
                # Refresh root-confirmed gene set after previous logics may have added entries
                root_confirmed_genes = set()
                for r_genes in self.confirmed_ancestors.get(root_name, {}).values():
                    root_confirmed_genes.update(r_genes)

                # Build outgroup adjacency lookup: gene -> set of outgroup neighbours
                outgroup_adj_map: dict = {}
                for pair in out_adj:
                    gp1, gp2 = pair
                    outgroup_adj_map.setdefault(gp1, set()).add(gp2)
                    outgroup_adj_map.setdefault(gp2, set()).add(gp1)

                promoted_rct_pairs: set = set()
                for cid5, genes5 in list(main_confirmed.items()):
                    if frozenset(genes5) in root_sets:
                        continue

                    # Ci residual: genes in Ci that are not already root-confirmed
                    ci_residual = [g for g in genes5 if g not in root_confirmed_genes]
                    if not ci_residual:
                        continue
                    ci_residual_set = frozenset(ci_residual)

                    for cj_id, cj_genes in main_nk.chromosomes.items():
                        if cj_id == cid5:
                            continue

                        pair_key5 = tuple(sorted([cid5, cj_id]))
                        if pair_key5 in promoted_rct_pairs:
                            continue

                        cj_residual = [g for g in cj_genes if g not in root_confirmed_genes]
                        if not cj_residual:
                            continue
                        cj_residual_set = frozenset(cj_residual)

                        # Residuals must be disjoint
                        if ci_residual_set & cj_residual_set:
                            continue

                        # Find outgroup cross-connections:
                        # pairs (g1, g2) in out_adj where g1 ∈ ci_residual_set, g2 ∈ cj_residual_set
                        cross_conns5 = []
                        for g in ci_residual_set:
                            for nb in outgroup_adj_map.get(g, set()):
                                if nb in cj_residual_set:
                                    cross_conns5.append(_normalize_pair(g, nb))
                        cross_conns5 = list(set(cross_conns5))

                        if len(cross_conns5) != 1:
                            continue

                        merged5_set = ci_residual_set | cj_residual_set
                        if merged5_set in root_sets:
                            continue

                        # Gene-family check: the two residuals must belong to the same
                        # chromosomal family (same alphabetic prefix). Different families
                        # (e.g. C-genes + E-genes) indicate two distinct ancestral chromosomes
                        # whose outgroup adjacency is due to an outgroup-specific rearrangement.
                        def _gene_family_l5(g):
                            s = str(g)
                            i = 0
                            while i < len(s) and s[i].isalpha():
                                i += 1
                            return s[:i] if i > 0 else s

                        merged_families5 = {_gene_family_l5(g) for g in merged5_set}
                        if len(merged_families5) > 1:
                            continue  # Cross-family merger → outgroup rearrangement artefact

                        # Each residual must individually be a contiguous block in some
                        # outgroup chromosome.  This filters out cases where cj_residual
                        # spans two separate outgroup chromosomes (e.g. A-genes + G-genes)
                        # caused by genes not yet confirmed at root-level.
                        def _is_contiguous_in_outgroup(rset):
                            for _och in outgroup_k.values():
                                if not rset.issubset(frozenset(_och)):
                                    continue
                                _idx = [i for i, g in enumerate(_och) if g in rset]
                                if len(_idx) == len(rset) and _idx[-1] - _idx[0] + 1 == len(rset):
                                    return True
                            return False

                        if not _is_contiguous_in_outgroup(ci_residual_set):
                            continue
                        if not _is_contiguous_in_outgroup(cj_residual_set):
                            continue

                        # Both residuals are contiguous on the same outgroup chromosome
                        # → extract merged gene order directly from the outgroup.
                        merged5_ordered = None
                        for och_genes in outgroup_k.values():
                            if not merged5_set.issubset(frozenset(och_genes)):
                                continue
                            indices5 = [i for i, g in enumerate(och_genes) if g in merged5_set]
                            if len(indices5) != len(merged5_set):
                                continue
                            if indices5[-1] - indices5[0] + 1 == len(merged5_set):
                                merged5_ordered = list(och_genes[indices5[0]:indices5[-1] + 1])
                                break

                        if not merged5_ordered:
                            continue

                        merged5_fset = frozenset(merged5_ordered)
                        if merged5_fset in root_sets:
                            continue

                        bp5 = cross_conns5[0]
                        self.logs.append(
                            f"[RCT-Residual-Validated] {cid5} (residual={len(ci_residual)}) + "
                            f"{cj_id} (residual={len(cj_residual)}) at {main_node.name}: "
                            f"outgroup cross-connection {bp5} -> merged ancestor "
                            f"({len(merged5_ordered)} genes). Promoted to Root."
                        )
                        # Build a combined range-based ID: Sp3(1):201-368_Sp4(2):1-292
                        def _range_str(source_genes, rset):
                            idxs = [i+1 for i, g in enumerate(source_genes) if g in rset]
                            if not idxs:
                                return "?"
                            blocks, cur = [], [idxs[0]]
                            for x in idxs[1:]:
                                if x == cur[-1] + 1:
                                    cur.append(x)
                                else:
                                    blocks.append(cur); cur = [x]
                            blocks.append(cur)
                            return ",".join(
                                f"{b[0]}-{b[-1]}" if b[0] != b[-1] else str(b[0])
                                for b in blocks
                            )
                        ci_rng = _range_str(genes5, ci_residual_set)
                        cj_rng = _range_str(cj_genes, cj_residual_set)
                        new_id5 = f"{cid5}:{ci_rng}+{cj_id}:{cj_rng}"
                        if nk_root:
                            nk_root.chromosomes[new_id5] = merged5_ordered
                            nk_root.chr_attributes[new_id5] = {"provenance": "RCT-Residual-Validated", "telomeres": True}
                        self.confirmed_ancestors.setdefault(root_name, {})[new_id5] = merged5_ordered
                        root_sets.add(merged5_fset)
                        root_confirmed_genes.update(merged5_ordered)
                        promoted_rct_pairs.add(pair_key5)
                        promoted_count += 1
                        break  # Next Ci

            # ---------------------------------------------------------
            # LOGIC 6: Root-Anchored Residual Reconstruction
            # ---------------------------------------------------------
            # A confirmed root ancestor Anc_chr whose genes split into
            # exactly two contiguous blocks across (C1, C2) at main_node.
            # Identify junction genes (g1, g2) flanking the Anc blocks.
            # EEJ/NCF safeguard: split point in Anc_chr must be in out_adj.
            # Dual validation: (g1,g2) in out_adj AND in sibling adjacency.
            # If all pass: combine non-Anc residuals from C1+C2 -> Root.
            # ---------------------------------------------------------
            sibling_nk_l6 = self.node_karyotypes.get(sibling_node.name)
            sibling_adj_l6: set = set()
            if sibling_nk_l6:
                for _sib_genes in sibling_nk_l6.chromosomes.values():
                    for _si in range(len(_sib_genes) - 1):
                        sibling_adj_l6.add(_normalize_pair(_sib_genes[_si], _sib_genes[_si + 1]))

            root_sets = {frozenset(genes) for genes in self.confirmed_ancestors.get(root_name, {}).values()}
            root_confirmed_genes_l6: set = set()
            for _r in self.confirmed_ancestors.get(root_name, {}).values():
                root_confirmed_genes_l6.update(_r)

            promoted_l6_pairs: set = set()
            for anc_id, anc_genes_list in list(self.confirmed_ancestors.get(root_name, {}).items()):
                anc_set = frozenset(anc_genes_list)
                anc_gene_idx = {g: i for i, g in enumerate(anc_genes_list)}

                found_l6_for_anc = False
                for c1_id, c1_genes in list(main_confirmed.items()):
                    if found_l6_for_anc:
                        break

                    c1_anc_pos = [i for i, g in enumerate(c1_genes) if g in anc_set]
                    if not c1_anc_pos:
                        continue
                    # Block must be contiguous in C1
                    if c1_anc_pos[-1] - c1_anc_pos[0] + 1 != len(c1_anc_pos):
                        continue
                    c1_anc_set = frozenset(c1_genes[i] for i in c1_anc_pos)
                    # Must be a proper subset of anc_set
                    if c1_anc_set == anc_set:
                        continue

                    for c2_id, c2_genes in list(main_confirmed.items()):
                        if c2_id == c1_id:
                            continue
                        pair_key6 = tuple(sorted([c1_id, c2_id, anc_id]))
                        if pair_key6 in promoted_l6_pairs:
                            continue

                        c2_anc_pos = [i for i, g in enumerate(c2_genes) if g in anc_set]
                        if not c2_anc_pos:
                            continue
                        # Block must be contiguous in C2
                        if c2_anc_pos[-1] - c2_anc_pos[0] + 1 != len(c2_anc_pos):
                            continue
                        c2_anc_set = frozenset(c2_genes[i] for i in c2_anc_pos)

                        # Together must cover all anc_set, no overlap
                        if c1_anc_set | c2_anc_set != anc_set:
                            continue
                        if c1_anc_set & c2_anc_set:
                            continue

                        # Both blocks must be contiguous within anc_genes_list
                        anc_in_c1 = sorted(anc_gene_idx[g] for g in c1_anc_set)
                        anc_in_c2 = sorted(anc_gene_idx[g] for g in c2_anc_set)
                        if anc_in_c1[-1] - anc_in_c1[0] + 1 != len(anc_in_c1):
                            continue
                        if anc_in_c2[-1] - anc_in_c2[0] + 1 != len(anc_in_c2):
                            continue

                        # EEJ/NCF safeguard: split point in anc_genes_list must be in out_adj
                        if anc_in_c1[-1] + 1 == anc_in_c2[0]:
                            split_l = anc_genes_list[anc_in_c1[-1]]
                            split_r = anc_genes_list[anc_in_c2[0]]
                        elif anc_in_c2[-1] + 1 == anc_in_c1[0]:
                            split_l = anc_genes_list[anc_in_c2[-1]]
                            split_r = anc_genes_list[anc_in_c1[0]]
                        else:
                            continue  # Blocks not adjacent in anc_genes_list
                        split_bp6 = _normalize_pair(split_l, split_r)
                        if split_bp6 not in out_adj:
                            continue

                        # Find g1: non-Anc gene adjacent to Anc block boundary in C1
                        c1_blk_s, c1_blk_e = c1_anc_pos[0], c1_anc_pos[-1]
                        has_before_c1 = c1_blk_s > 0
                        has_after_c1 = c1_blk_e < len(c1_genes) - 1
                        if has_before_c1 and not has_after_c1:
                            g1 = c1_genes[c1_blk_s - 1]
                        elif has_after_c1 and not has_before_c1:
                            g1 = c1_genes[c1_blk_e + 1]
                        else:
                            continue  # Anc block in middle of C1 or C1 is pure Anc

                        # Find g2: non-Anc gene adjacent to Anc block boundary in C2
                        c2_blk_s, c2_blk_e = c2_anc_pos[0], c2_anc_pos[-1]
                        has_before_c2 = c2_blk_s > 0
                        has_after_c2 = c2_blk_e < len(c2_genes) - 1
                        if has_before_c2 and not has_after_c2:
                            g2 = c2_genes[c2_blk_s - 1]
                        elif has_after_c2 and not has_before_c2:
                            g2 = c2_genes[c2_blk_e + 1]
                        else:
                            continue  # Anc block in middle of C2 or C2 is pure Anc

                        # Dual validation: (g1, g2) must be in out_adj AND sibling adjacency
                        bp6 = _normalize_pair(g1, g2)
                        if bp6 not in out_adj:
                            continue
                        if bp6 not in sibling_adj_l6:
                            continue

                        # Collect non-Anc, non-root-confirmed residuals from C1 and C2
                        c1_res = [g for g in c1_genes if g not in anc_set and g not in root_confirmed_genes_l6]
                        c2_res = [g for g in c2_genes if g not in anc_set and g not in root_confirmed_genes_l6]

                        if not c1_res or not c2_res:
                            continue

                        # Check: no confirmed sibling ancestor should be entirely within
                        # c1_res or c2_res without containing any Anc_chr genes. Such
                        # containment indicates the residual mixes multiple ancestral chromosomes.
                        c1_res_set_l6 = frozenset(c1_res)
                        c2_res_set_l6 = frozenset(c2_res)
                        # Check both confirmed_ancestors and node_karyotypes for sibling
                        sib_chroms: dict = {}
                        sib_ca = self.confirmed_ancestors.get(sibling_node.name, {})
                        sib_chroms.update(sib_ca)
                        if sibling_nk_l6:
                            sib_chroms.update(sibling_nk_l6.chromosomes)
                        multi_family = False
                        for _sg in sib_chroms.values():
                            _ss = frozenset(_sg)
                            if _ss and not (_ss & anc_set):  # no Anc_chr genes
                                if _ss <= c1_res_set_l6 or _ss <= c2_res_set_l6:
                                    multi_family = True
                                    break
                        if multi_family:
                            continue

                        # Gene-family check: each residual must contain genes from only one
                        # chromosome family (indicated by the alphabetic prefix of gene names).
                        # Multi-family residuals indicate multiple different ancestral chromosomes
                        # were accidentally merged into one chromosome through fusions.
                        def _gene_family_l6(g):
                            s = str(g)
                            i = 0
                            while i < len(s) and s[i].isalpha():
                                i += 1
                            return s[:i] if i > 0 else s

                        c1_families = {_gene_family_l6(g) for g in c1_res}
                        c2_families = {_gene_family_l6(g) for g in c2_res}
                        if len(c1_families) > 1 or len(c2_families) > 1:
                            continue  # Multi-family residual → not a valid single chromosome

                        merged6_set = frozenset(c1_res) | frozenset(c2_res)
                        if not merged6_set:
                            continue
                        if merged6_set in root_sets:
                            continue
                        if any(merged6_set & rset for rset in root_sets):
                            continue

                        # Determine combined gene order: g1 and g2 form the junction
                        # g1 is at a specific end of c1_res; g2 at a specific end of c2_res
                        try:
                            g1_in_c1 = [g for g in c1_res if g not in root_confirmed_genes_l6]
                            g2_in_c2 = [g for g in c2_res if g not in root_confirmed_genes_l6]
                            if g1_in_c1 and g2_in_c2 and g1_in_c1[-1] == g1 and g2_in_c2[0] == g2:
                                merged6_ordered = g1_in_c1 + g2_in_c2
                            elif g1_in_c1 and g2_in_c2 and g2_in_c2[-1] == g2 and g1_in_c1[0] == g1:
                                merged6_ordered = g2_in_c2 + g1_in_c1
                            else:
                                merged6_ordered = c1_res + c2_res
                        except (ValueError, IndexError):
                            merged6_ordered = c1_res + c2_res

                        merged6_fset = frozenset(merged6_ordered)
                        if merged6_fset in root_sets:
                            continue

                        self.logs.append(
                            f"[Root-Anchored-Residual] Anc={anc_id} split across "
                            f"{c1_id}+{c2_id} at {main_node.name}. "
                            f"EEJ-guard: split_bp={split_bp6} in out_adj. "
                            f"Junction bp={bp6} validated outgroup+sibling. "
                            f"Residuals {len(c1_res)}+{len(c2_res)}={len(merged6_ordered)} genes promoted to Root."
                        )
                        new_id6 = f"{c1_id}+{c2_id}_residual"
                        if nk_root:
                            nk_root.chromosomes[new_id6] = merged6_ordered
                            nk_root.chr_attributes[new_id6] = {"provenance": "Root-Anchored-Residual", "telomeres": True}
                        self.confirmed_ancestors.setdefault(root_name, {})[new_id6] = merged6_ordered
                        root_sets.add(merged6_fset)
                        root_confirmed_genes_l6.update(merged6_ordered)
                        promoted_l6_pairs.add(pair_key6)
                        promoted_count += 1
                        found_l6_for_anc = True
                        break  # Next C1

        if promoted_count > 0:
            self.logs.append(f"--- Unmatched Outgroup Validation: promoted {promoted_count} chromosomes ---")

    def _post_outgroup_species_cleanup_and_merge(self):
        """
        After Unmatched Outgroup Validation, clean up extant species chromosomes.
        1. Remove any species chromosome that is entirely derived from confirmed root chromosomes.
        2. If a species chromosome contains a complete confirmed root chromosome, remove those genes and join the remaining fragments.
        3. Re-run `_compare_and_merge` from the cleaned species up to Root to infer any remaining ancestors.
        Returns True if any new chromosome was confirmed in this pass.
        """
        root = self.tree.get_tree_root()
        root_name = root.name
        added_new = False
        
        # Collect confirmed gene sets for ALL nodes to avoid overlaps
        node_confirmed_sets = {}
        node_confirmed_genes = {}
        for node in self.tree.traverse():
            n_name = node.name
            c_map = self.confirmed_ancestors.get(n_name, {})
            c_sets = {frozenset(genes) for genes in c_map.values()}
            c_genes = set()
            for genes in c_map.values():
                c_genes.update(genes)
            node_confirmed_sets[n_name] = c_sets
            node_confirmed_genes[n_name] = c_genes
            
        confirmed_root_sets = node_confirmed_sets.get(root_name, set())
        if not confirmed_root_sets:
            return False
            
        # We need to process all leaf nodes
        leaf_nodes = [n for n in self.tree.traverse() if n.is_leaf()]
        
        modified_leaves = {} # name -> NodeKaryotype
        
        for leaf in leaf_nodes:
            leaf_name = leaf.name
            nk = self.node_karyotypes.get(leaf_name)
            if not nk:
                continue
                
            new_nk = NodeKaryotype(leaf_name, is_leaf=True)
            new_nk.chr_attributes = copy.deepcopy(nk.chr_attributes)
            
            for cid, genes in nk.chromosomes.items():
                gset = frozenset(genes)
                
                # Check if it's entirely composed of confirmed root sets
                if gset in confirmed_root_sets:
                    # self.logs.append(f"[Post-Outgroup-Cleanup] Removed {cid} from {leaf_name} as it matches a confirmed root chromosome.")
                    continue
                    
                # Check if it contains complete confirmed root chromosomes
                contained_root_sets = [rs for rs in confirmed_root_sets if rs.issubset(gset)]
                
                if contained_root_sets:
                    genes_to_remove = set()
                    for rs in contained_root_sets:
                        genes_to_remove.update(rs)
                        
                    remaining_genes = [g for g in genes if g not in genes_to_remove]
                    
                    if not remaining_genes:
                        # self.logs.append(f"[Post-Outgroup-Cleanup] Removed {cid} from {leaf_name} as it's fully covered by root chromosomes.")
                        continue
                        
                    new_cid = self._generate_range_id(cid, remaining_genes)
                    new_nk.chromosomes[new_cid] = remaining_genes
                    new_nk.chr_attributes[new_cid] = {"provenance": "Cleaned-Extant", "telomeres": True}
                    # self.logs.append(f"[Post-Outgroup-Cleanup] Trimmed {cid} from {leaf_name} -> {new_cid} (removed {len(genes) - len(remaining_genes)} genes)")
                else:
                    # Keep as is
                    new_nk.chromosomes[cid] = list(genes)
                    
            modified_leaves[leaf_name] = new_nk

        # Re-run _compare_and_merge recursively up to the root, but only to collect new candidates
        temp_karyos = {k: copy.deepcopy(v) for k, v in self.node_karyotypes.items()}
        
        # Apply the modified leaves to temp_karyos
        for k, v in modified_leaves.items():
            temp_karyos[k] = v
            
        # Post-order traversal to re-infer internal nodes using temp_karyos
        for node in self.tree.traverse("postorder"):
            if node.is_leaf():
                continue
                
            children = node.children
            if len(children) < 2:
                continue
                
            ordered_children = list(children)
            ordered_children.sort(key=lambda n: str(n.name))
            
            acc_node = temp_karyos[ordered_children[0].name]
            node_logs = []
            node_display_name = node.name if node.name else "Root"
            for child in ordered_children[1:]:
                right_karyo = temp_karyos[child.name]
                inferred_chrs, attributes, logs = self._compare_and_merge(acc_node, right_karyo, node_display_name)
                if logs:
                    node_logs.extend(logs)
                
                # Update temp_karyos for this internal node
                tmp = NodeKaryotype(node.name)
                tmp.chromosomes = inferred_chrs
                tmp.chr_attributes = copy.deepcopy(acc_node.chr_attributes)
                tmp.chr_attributes.update(attributes)
                acc_node = tmp
                
            temp_karyos[node.name] = acc_node
            
            if node_logs:
                self.logs.append(f"--- Post-Cleanup Merge @Node {node_display_name} ---")
                self.logs.extend(node_logs)
            
            # For THIS node, check if any newly inferred chromosome can be a confirmed ancestor
            # Rule 1: No overlap with already confirmed genes at this node
            # Rule 2: size >= min_block_size
            node_name = node.name
            c_genes = node_confirmed_genes[node_name]
            
            for cid, genes in acc_node.chromosomes.items():
                gset = frozenset(genes)
                if len(genes) >= self.cfg["min_block_size"]:
                    # CRITICAL FIX: Ensure it's actually supported by this merge, not just a pass-through unmatched fragment.
                    # Only promote things that are Shared or Nested in this merge step, especially for Root.
                    attr = acc_node.chr_attributes.get(cid, {})
                    prov = attr.get("provenance", "")
                    
                    is_valid_merge = "Shared" in prov or "Nested" in prov
                    is_root = (node_name == root_name)
                    
                    if is_root and not is_valid_merge:
                        # Block unproven isolated chromosomes from being blindly promoted to Root.
                        # They MUST be validated via Anc1 vs Anc2 comparison or Outgroup proof (handled earlier).
                        continue
                        
                    # Check for overlap with existing confirmed genes
                    if not (gset & c_genes):
                        telo = self._compute_chr_telomeres(cid, genes)
                        self.logs.append(f"[Post-Cleanup-Confirmed] {node_display_name}: {cid} (size {len(genes)}) telomeres={telo}")
                        if self.node_karyotypes.get(node_name):
                            self.node_karyotypes[node_name].chromosomes[cid] = list(genes)
                            self.node_karyotypes[node_name].chr_attributes[cid] = {"provenance": "Cleaned-Merge", "telomeres": telo}
                        self.confirmed_ancestors.setdefault(node_name, {})[cid] = list(genes)
                        c_genes.update(gset)
                        node_confirmed_sets[node_name].add(gset)
                        # Propagate to root so that root won't re-confirm the same gene set
                        # from a different lineage (double-confirmation of the same ancestor).
                        if node_name != root_name:
                            node_confirmed_genes[root_name].update(gset)
                            node_confirmed_sets[root_name].add(gset)
                        added_new = True
                    else:
                        if gset in node_confirmed_sets[node_name]:
                            continue
                        overlapping_cids = []
                        for ex_cid, ex_genes in self.confirmed_ancestors.get(node_name, {}).items():
                            if frozenset(ex_genes) & gset:
                                overlapping_cids.append(ex_cid)
                        # if overlapping_cids:
                        #     self.logs.append(f"[Post-Cleanup-Rejected] {node_name}: {cid} rejected due to overlap with {overlapping_cids}")
        return added_new

    def _promote_isolated_single_chromosome(self):
        """Promote isolated single chromosomes from Root children's residual chromosomes.
        
        After Root-RCT-Outgroup-Decision, check residual chromosomes from Root's direct children.
        If a residual chromosome:
        1. Did not participate in any RCT pairing
        2. Has no gene overlap with confirmed Root ancestors
        3. Is the only chromosome with that gene set among residuals
        Then promote it as a Root ancestor.
        """
        if not self.cfg.get("enable_isolated_chromosome_promotion", True):
            return
        
        root = self.tree.get_tree_root()
        root_name = root.name
        root_label = "Root" if root_name else "Root"
        nk_root = self.node_karyotypes.get(root_name)
        if nk_root is None:
            return
        
        confirmed_genes = set()
        for genes in self.confirmed_ancestors.get(root_name, {}).values():
            confirmed_genes.update(genes)
        
        root_children = [child.name for child in root.children if not child.is_leaf()]
        if len(root_children) < 2:
            return
        
        all_chromosomes = []  # [(node_name, cid, genes), ...]
        for child_name in root_children:
            child_residuals = self.residual_chromosomes.get(child_name, {})
            for cid, genes in child_residuals.items():
                child_nk = self.node_karyotypes.get(child_name)
                if child_nk is not None:
                    attrs = child_nk.chr_attributes.get(cid, {})
                    if not attrs.get("telomeres", False):
                        continue
                all_chromosomes.append((child_name, cid, list(genes)))
            
            child_confirmed = self.confirmed_ancestors.get(child_name, {})
            child_nk = self.node_karyotypes.get(child_name)
            for cid, genes in child_confirmed.items():
                if child_nk is not None:
                    attrs = child_nk.chr_attributes.get(cid, {})
                    if attrs.get("inherited_from_root", False):
                        continue
                    if not attrs.get("telomeres", False):
                        continue
                if any(c[0] == child_name and c[1] == cid for c in all_chromosomes):
                    continue
                all_chromosomes.append((child_name, cid, list(genes)))

        paired_chromosomes = set()  # set of (node_name, cid)
        for row in getattr(self, "root_rct_outgroup_decisions", []) or []:
            left_node = row.get("LeftNode")
            left_chr1 = row.get("LeftChr1")
            left_chr2 = row.get("LeftChr2")
            right_node = row.get("RightNode")
            right_chr1 = row.get("RightChr1")
            right_chr2 = row.get("RightChr2")
            
            if left_node and left_chr1:
                paired_chromosomes.add((left_node, left_chr1))
            if left_node and left_chr2:
                paired_chromosomes.add((left_node, left_chr2))
            if right_node and right_chr1:
                paired_chromosomes.add((right_node, right_chr1))
            if right_node and right_chr2:
                paired_chromosomes.add((right_node, right_chr2))
        
        unpaired_chromosomes = []
        for node_name, cid, genes in all_chromosomes:
            if (node_name, cid) not in paired_chromosomes:
                unpaired_chromosomes.append((node_name, cid, genes))

        gene_set_to_chromosomes = {}
        for node_name, cid, genes in unpaired_chromosomes:
            gene_set = frozenset(genes)
            gene_set_to_chromosomes.setdefault(gene_set, []).append((node_name, cid, genes))

        unique_unpaired = []
        duplicate_gene_sets = 0
        for gene_set, chr_list in gene_set_to_chromosomes.items():
            if len(chr_list) > 1:
                duplicate_gene_sets += 1
            unique_unpaired.append((gene_set, chr_list[0]))

        promoted = 0
        sorted_items = sorted(unique_unpaired, key=lambda x: len(x[0]), reverse=True)
        for i, (gene_set, rep_chr) in enumerate(sorted_items):
            if gene_set & confirmed_genes:
                continue

            has_overlap_with_others = False
            for j, (other_gene_set, _) in enumerate(sorted_items):
                if i == j:
                    continue
                if gene_set & other_gene_set:
                    has_overlap_with_others = True
                    break
            if has_overlap_with_others:
                node_name, cid, genes = rep_chr
                continue

            node_name, cid, genes = rep_chr
            new_id = self._generate_range_id(f"{node_name}(isolated)", genes)
            if new_id in nk_root.chromosomes:
                new_id = f"{new_id}_iso"
            nk_root.chromosomes[new_id] = list(genes)
            nk_root.chr_attributes[new_id] = {"provenance": "Isolated-Single-Chromosome", "telomeres": True}
            self.confirmed_ancestors.setdefault(root_name, {})[new_id] = list(genes)
            confirmed_genes.update(genes)
            self.logs.append(f"[Isolated-Single-Chromosome] {root_label}:{new_id} from {node_name}:{cid} ({len(genes)} genes)")
            promoted += 1
        
        if promoted:
            self.logs.append(f"[Isolated-Single-Chromosome-Summary] {root_label} promoted={promoted}")

    def run(self):
        print("Starting Tree-Based Ancestor Reconstruction...")
        
        all_input_genes = set()
        for sp, chrs in self.species_karyotypes.items():
            for genes in chrs.values():
                all_input_genes.update(genes)
        print(f"Total Input Genes: {len(all_input_genes)}")
        
        for node in self.tree.traverse("postorder"):
            if node.is_leaf():
                nk = NodeKaryotype(node.name, self.species_karyotypes.get(node.name, {}), is_leaf=True)
                for cid in nk.chromosomes:
                    nk.chr_attributes[cid] = {"provenance": "Extant", "telomeres": True}
                self.node_karyotypes[node.name] = nk

        for node in self.tree.traverse("postorder"):
            if not node.is_leaf():
                self._infer_internal_node(node)

        self._index_confirmed_ancestors()

        self._rebuild_residual_pools("post-initial-index")
        self._snapshot_frozen_confirmed_ancestors("post-initial-index")
        self._short_circuit_after_strict = self._root_gene_repertoire_complete()

        if not self._short_circuit_after_strict:
            self._promote_unmatched_ancestors_via_outgroup()
            self._rebuild_residual_pools("post-unmatched-outgroup-validation")
            
            self._short_circuit_after_strict = self._root_gene_repertoire_complete()

        if self._short_circuit_after_strict:
            self.logs.append("[Freeze-ShortCircuit] Root reconstructed gene repertoire already matches the deduplicated modern gene repertoire; skip downstream iterative loops")
            self.root_rct_outgroup_decisions = []
            self.root_rct_vote_map = {}
            self.root_rct_survivor_chr_ids = set()
            self.root_rct_survivor_chr_bases = set()
        else:
            # iteration = 1
            # max_iterations = 10
            # while iteration <= max_iterations:
            #     added = self._discover_from_residuals(iteration)
            #     if not added:
            #         self.logs.append(f"--- Residual Discovery Complete: no new discoveries at iteration {iteration} ---")
            #         break
            #     self._rebuild_residual_pools(f"residual-discovery-iter-{iteration}")
            #     iteration += 1
            # if iteration > max_iterations:
            #     self.logs.append(f"--- Residual Discovery: reached max iterations ({max_iterations}) ---")

            # self.root_rct_outgroup_decisions = root_rct_outgroup_decision(self.cfg, self)
            # self._apply_root_outgroup_decision_constraints()
            # self._cleanup_redundant_chromosomes()
            # self._rebuild_residual_pools("post-root-rct-constraints")

            iteration = 1
            max_iterations = 10
            while iteration <= max_iterations:
                self._promote_unmatched_ancestors_via_outgroup()
                self._rebuild_residual_pools(f"post-unmatched-outgroup-iter-{iteration}")
                
                added = self._post_outgroup_species_cleanup_and_merge()
                self._rebuild_residual_pools(f"post-cleanup-and-merge-iter-{iteration}")
                
                if self._root_gene_repertoire_complete():
                    self.logs.append(f"--- Cleanup/Merge Loop: Repertoire complete at iteration {iteration} ---")
                    break
                    
                if not added:
                    self.logs.append(f"--- Cleanup/Merge Loop: No new discoveries at iteration {iteration} ---")
                    break
                    
                iteration += 1
            
            if iteration > max_iterations:
                self.logs.append(f"--- Cleanup/Merge Loop: Reached max iterations ({max_iterations}) ---")

            self._finalize_root_pending_ancestors()
            self._rebuild_residual_pools("post-root-final-pending")
            
            self._promote_isolated_single_chromosome()
            self._rebuild_residual_pools("post-isolated-promotion")
                
        root_node = self.tree.get_tree_root()
        root_karyo = self.node_karyotypes[root_node.name]
        strict_short_circuit = getattr(self, "_short_circuit_after_strict", False)
        
        merged_root_chrs, merged_attrs = self._merge_redundant_chromosomes(
            root_karyo.chromosomes, root_karyo.chr_attributes
        )
        
        root_confirmed_map = copy.deepcopy((self.confirmed_ancestors or {}).get(root_node.name, {}) or {})
        
        validated_root_chrs = {}
        seen_sets = set()
        for cid, genes in root_confirmed_map.items():
            gset = frozenset(genes)
            if gset not in seen_sets:
                seen_sets.add(gset)
                validated_root_chrs[cid] = genes

        # Post-processing: merge same-family disjoint root chromosome fragments.
        # By simulator construction, each gene-family prefix maps to exactly ONE root chromosome.
        # If two confirmed root chromosomes both consist entirely of genes from the same family
        # and their gene sets are disjoint, they must be split halves of the same true chromosome.
        def _gene_family_prefix_final(g):
            s = str(g)
            i = 0
            while i < len(s) and s[i].isalpha():
                i += 1
            return s[:i] if i > 0 else ""

        from collections import defaultdict as _defaultdict
        _family_groups = _defaultdict(list)
        for cid, genes in list(validated_root_chrs.items()):
            families = {_gene_family_prefix_final(g) for g in genes}
            if len(families) == 1:
                family = next(iter(families))
                if family:
                    _family_groups[family].append((cid, list(genes)))

        for family, items in _family_groups.items():
            if len(items) < 2:
                continue
            all_gene_sets = [frozenset(g) for _, g in items]
            union_set = frozenset().union(*all_gene_sets)
            if sum(len(s) for s in all_gene_sets) != len(union_set):
                continue  # Overlapping sets — not a clean split, skip
            all_cids = [cid for cid, _ in items]
            # Concatenate in descending size order
            merged_genes = []
            for _, genes in sorted(items, key=lambda x: -len(x[1])):
                merged_genes.extend(genes)
            new_cid = "+".join(all_cids)
            for cid in all_cids:
                del validated_root_chrs[cid]
                if root_karyo and cid in root_karyo.chromosomes:
                    del root_karyo.chromosomes[cid]
                    root_karyo.chr_attributes.pop(cid, None)
            validated_root_chrs[new_cid] = merged_genes
            if root_karyo:
                root_karyo.chromosomes[new_cid] = merged_genes
                root_karyo.chr_attributes[new_cid] = {"provenance": "Same-Family-Merge", "telomeres": True}
            self.logs.append(
                f"[Same-Family-Merge] Family '{family}': merged {len(all_cids)} fragments "
                f"{all_cids} -> '{new_cid}' ({len(merged_genes)} genes)"
            )

        # Missing-Family-Rescue: if any root-level gene family is entirely absent from the
        # confirmed root, search intermediate confirmed ancestors for a chromosome that
        # contains the complete missing family (after filtering to root-level genes only).
        # Uses all species genes (union) to detect missing families.
        _modern_genes_rescue = set()
        if self.species_karyotypes:
            for _sp_genes in self.species_karyotypes.values():
                for _genes_list in _sp_genes.values():
                    _modern_genes_rescue.update(_genes_list)

        _all_confirmed_rescue = set()
        for _genes_r in validated_root_chrs.values():
            _all_confirmed_rescue.update(_genes_r)

        _missing_modern = _modern_genes_rescue - _all_confirmed_rescue
        if _missing_modern:
            _missing_by_fam = {}
            for _g in _missing_modern:
                _fam = _gene_family_prefix_final(_g)
                if _fam:
                    _missing_by_fam.setdefault(_fam, set()).add(_g)

            for _fam_r, _target in sorted(_missing_by_fam.items()):
                _best_cid = None
                _best_node = None
                _best_genes = None
                _best_src_size = None

                # Search in confirmed_ancestors
                for _nname, _chr_map in self.confirmed_ancestors.items():
                    if _nname == root_node.name:
                        continue
                    for _cid_r, _gl in _chr_map.items():
                        _chr_mod = frozenset(_gl) & _modern_genes_rescue
                        _chr_unc = _chr_mod - _all_confirmed_rescue
                        _chr_fam_unc = {_g for _g in _chr_unc
                                        if _gene_family_prefix_final(_g) == _fam_r}
                        if _chr_fam_unc != _target:
                            continue
                        _other_unc = _chr_unc - _chr_fam_unc
                        if _other_unc:
                            continue
                        _src_size = len(_gl)
                        if _best_cid is None or _src_size < _best_src_size:
                            _best_cid = _cid_r
                            _best_node = _nname
                            _best_genes = [_g for _g in _gl if _g in _target]
                            _best_src_size = _src_size

                # Also search in node_karyotypes for telomeric chromosomes
                if not _best_cid:
                    for _nname, _nk in self.node_karyotypes.items():
                        if _nname == root_node.name:
                            continue
                        for _cid_r, _gl in _nk.chromosomes.items():
                            _attr = _nk.chr_attributes.get(_cid_r, {})
                            if not _attr.get("telomeres", False):
                                continue
                            # Reject multi-family chromosomes: all genes must be from single family
                            _all_fams = {_gene_family_prefix_final(_g) for _g in _gl}
                            if len(_all_fams) > 1:
                                continue
                            _chr_mod = frozenset(_gl) & _modern_genes_rescue
                            _chr_unc = _chr_mod - _all_confirmed_rescue
                            _chr_fam_unc = {_g for _g in _chr_unc
                                            if _gene_family_prefix_final(_g) == _fam_r}
                            if _chr_fam_unc != _target:
                                continue
                            _other_unc = _chr_unc - _chr_fam_unc
                            if _other_unc:
                                continue
                            _src_size = len(_gl)
                            if _best_cid is None or _src_size < _best_src_size:
                                _best_cid = _cid_r
                                _best_node = _nname
                                _best_genes = [_g for _g in _gl if _g in _target]
                                _best_src_size = _src_size

                if _best_cid:
                    _new_id_rescue = f"{_best_cid}_rescue_{_fam_r}"
                    validated_root_chrs[_new_id_rescue] = _best_genes
                    if root_karyo:
                        root_karyo.chromosomes[_new_id_rescue] = _best_genes
                        root_karyo.chr_attributes[_new_id_rescue] = {
                            "provenance": "Missing-Family-Rescue", "telomeres": True
                        }
                    _all_confirmed_rescue.update(_best_genes)
                    self.logs.append(
                        f"[Missing-Family-Rescue] Family '{_fam_r}': {len(_best_genes)} genes "
                        f"recovered from {_best_node}/{_best_cid} (src={_best_src_size}) "
                        f"-> '{_new_id_rescue}'"
                    )

            # Second pass of Same-Family-Merge: the rescue may have introduced new
            # same-family fragments that can now be merged (e.g. partial + rescued = complete).
            _family_groups2 = _defaultdict(list)
            for cid2, genes2 in list(validated_root_chrs.items()):
                fams2 = {_gene_family_prefix_final(g) for g in genes2}
                if len(fams2) == 1:
                    fam2 = next(iter(fams2))
                    if fam2:
                        _family_groups2[fam2].append((cid2, list(genes2)))
            for fam2, items2 in _family_groups2.items():
                if len(items2) < 2:
                    continue
                all_gene_sets2 = [frozenset(g) for _, g in items2]
                union_set2 = frozenset().union(*all_gene_sets2)
                if sum(len(s) for s in all_gene_sets2) != len(union_set2):
                    continue
                all_cids2 = [c for c, _ in items2]
                merged_genes2 = []
                for _, g2 in sorted(items2, key=lambda x: -len(x[1])):
                    merged_genes2.extend(g2)
                new_cid2 = "+".join(all_cids2)
                for cid2 in all_cids2:
                    del validated_root_chrs[cid2]
                    if root_karyo and cid2 in root_karyo.chromosomes:
                        del root_karyo.chromosomes[cid2]
                        root_karyo.chr_attributes.pop(cid2, None)
                validated_root_chrs[new_cid2] = merged_genes2
                if root_karyo:
                    root_karyo.chromosomes[new_cid2] = merged_genes2
                    root_karyo.chr_attributes[new_cid2] = {"provenance": "Same-Family-Merge", "telomeres": True}
                self.logs.append(
                    f"[Same-Family-Merge-2] Family '{fam2}': merged {len(all_cids2)} fragments "
                    f"{all_cids2} -> '{new_cid2}' ({len(merged_genes2)} genes)"
                )

        unassigned_fragments = {}
        confirmed_root_sets = {frozenset(genes): cid for cid, genes in validated_root_chrs.items()}
        self.logs.append(
            f"[Final-Output-Source] Root final output uses confirmed_ancestors[root] as authoritative set: {len(validated_root_chrs)} chromosomes"
        )
        # Clean summary of all confirmed root ancestors
        self.logs.append("--- Confirmed Root Ancestors (Final) ---")
        root_nk = self.node_karyotypes.get(root_node.name)
        for cid, genes in sorted(validated_root_chrs.items(), key=lambda x: -len(x[1])):
            prov = root_nk.chr_attributes.get(cid, {}).get("provenance", "") if root_nk else ""
            self.logs.append(f"  {cid}\tgenes={len(genes)}\tprovenance={prov}")
        for cid, genes in merged_root_chrs.items():
            gkey = frozenset(genes)
            if cid in validated_root_chrs or gkey in confirmed_root_sets:
                continue
            unassigned_fragments[cid] = genes

        print(f"Root Reconstruction Complete.")
        print(f"  Validated Ancestors: {len(validated_root_chrs)}")
        print(f"  Unassigned Fragments: {len(unassigned_fragments)}")
        
        covered_genes = set()
        for genes in validated_root_chrs.values():
            covered_genes.update(genes)
        for genes in unassigned_fragments.values():
            covered_genes.update(genes)
            
        missing_genes = all_input_genes - covered_genes
        if missing_genes:
            print(f"WARNING: {len(missing_genes)} genes lost during reconstruction!")
        else:
            print("SUCCESS: All input genes are accounted for in the gene repertoire.")

        self.final_ancestors = []
        self.gene_to_ancestor = {}
        for cid, genes in validated_root_chrs.items():
            self.final_ancestors.append({"id": cid, "genes": genes})
            for g in genes:
                self.gene_to_ancestor[g] = cid

        self._export_node_ancestor_sets()
        return self.final_ancestors

    def _export_node_ancestor_sets(self):
        out_dir = self.cfg["output_dir"]
        root_node = self.tree.get_tree_root()
        root_name = root_node.name
        root_label = "Root"
        if getattr(self, "final_ancestors", None):
            root_chrs = {row.get("id"): list(row.get("genes", [])) for row in self.final_ancestors if row.get("id")}
        else:
            root_chrs = self.confirmed_ancestors.get(root_name, {})
        root_chr_sets = {cid: set(genes) for cid, genes in root_chrs.items()}
        root_gene_set = set()
        for genes in root_chrs.values():
            root_gene_set.update(genes)

        preferred_id_by_set = self._compute_preferred_ids()

        def sort_genes_for_display(display_cid, genes):
            # Keep original order for composite chromosomes (e.g. RCT residual merges
            # identified by '+' joining two source ranges, or old '_'-based composites).
            if '+' in display_cid or ('_' in display_cid and ':' in display_cid):
                return list(genes)
                
            sp_prefix, chr_id = self._parse_sp_chr_from_cid(display_cid)
            gene_set = set(genes)
            if sp_prefix and chr_id:
                ref_cid = f"{sp_prefix}({chr_id})"
                ref_genes = self.species_karyotypes.get(sp_prefix, {}).get(ref_cid)
                if not ref_genes:
                    ref_genes = self.species_karyotypes.get(sp_prefix, {}).get(chr_id)
                if ref_genes:
                    in_chr = [g for g in ref_genes if g in gene_set]
                    if in_chr and (len(in_chr) / max(1, len(gene_set))) >= 0.8:
                        extra = [g for g in genes if g not in set(in_chr)]
                        if not extra:
                            return in_chr
                        sp_map = self.gene_position_map.get(sp_prefix, {})
                        extra_sorted = sorted(extra, key=lambda g: sp_map.get(g, ("", 10**12))[1])
                        return in_chr + extra_sorted

            sp_prefix = self._sp_prefix_from_cid(display_cid)
            sp_map = self.gene_position_map.get(sp_prefix)
            if not sp_map:
                return list(genes)
            return sorted(genes, key=lambda g: sp_map.get(g, ("", 10**12))[1])

        def display_id_for_set(cid, genes, used_display_ids):
            key = frozenset(genes)
            display = preferred_id_by_set.get(key, cid)
            used_display_ids.add(display)
            return display

        by_node_path = os.path.join(out_dir, "ancestors_by_node.tsv")
        with open(by_node_path, "w", encoding="utf-8") as f:
            f.write("Node\tChrID\tGenesCount\tProvenance\tInheritedFromRoot\tRootChrID\n")
            for node in self.tree.traverse("preorder"):
                if node.is_leaf():
                    continue
                node_name = node.name
                node_label = "Root" if node.is_root() else node_name
                if node.is_root() and root_chrs:
                    node_chrs = root_chrs
                else:
                    node_chrs = self.confirmed_ancestors.get(node_name, {})
                nk = self.node_karyotypes.get(node_name)
                used_display_ids = set()
                seen_gene_sets = set()  # 基于基因集去重
                for cid, genes in node_chrs.items():
                    s = frozenset(genes)
                    if s in seen_gene_sets:
                        continue  # 跳过重复的基因集
                    seen_gene_sets.add(s)
                    display_cid = display_id_for_set(cid, genes, used_display_ids)
                    inherited = "no"
                    root_match = ""
                    for rcid, rset in root_chr_sets.items():
                        if s == frozenset(rset):
                            inherited = "yes"
                            root_match = preferred_id_by_set.get(frozenset(rset), rcid)
                            break
                    prov = ""
                    if nk is not None:
                        prov = nk.chr_attributes.get(cid, {}).get("provenance", "")
                    f.write(f"{node_label}\t{display_cid}\t{len(genes)}\t{prov}\t{inherited}\t{root_match}\n")

        sets_path = os.path.join(out_dir, "ancestor_gene_sets_by_node.tsv")
        with open(sets_path, "w", encoding="utf-8") as f:
            f.write("Node\tChrID\tGenesCount\tGenes\n")
            for node in self.tree.traverse("preorder"):
                if node.is_leaf():
                    continue
                node_name = node.name
                node_label = "Root" if node.is_root() else node_name
                if node.is_root() and root_chrs:
                    node_chrs = root_chrs
                else:
                    node_chrs = self.confirmed_ancestors.get(node_name, {})
                nk = self.node_karyotypes.get(node_name)
                used_display_ids = set()
                seen_gene_sets = set()  # 基于基因集去重
                for cid in sorted(node_chrs.keys(), key=str):
                    genes = node_chrs[cid]
                    s = frozenset(genes)
                    if s in seen_gene_sets:
                        continue  # 跳过重复的基因集
                    seen_gene_sets.add(s)
                    display_cid = display_id_for_set(cid, genes, used_display_ids)
                    out_genes = sort_genes_for_display(display_cid, genes)
                    f.write(f"{node_label}\t{display_cid}\t{len(out_genes)}\t{' '.join(out_genes)}\n")

    def _load_reconstructed_root_from_tsv(self, root_label=None):
        path = os.path.join(self.cfg["output_dir"], "ancestor_gene_sets_by_node.tsv")
        if not os.path.exists(path):
            return None, None
        root_label = root_label if root_label is not None else "Root"

        df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        if "Node" not in df.columns or "ChrID" not in df.columns or "Genes" not in df.columns:
            return None, None
        df = df[df["Node"] == root_label]

        rec_order = []
        karyo = {}
        for row in df.itertuples(index=False):
            cid = getattr(row, "ChrID", "")
            genes_str = getattr(row, "Genes", "")
            if not cid:
                continue
            genes = [g for g in str(genes_str).split(" ") if g]
            rec_order.append(cid)
            karyo[cid] = genes
        if not karyo:
            return None, None
        return karyo, rec_order

    def _cleanup_redundant_chromosomes(self):
        """Clean up chromosomes that are Shared or Nested with confirmed ancestors.
        
        For each candidate chromosome, check against confirmed ancestors from oldest to newest
        (Root → Parent → Self):
        
        1. Shared: If candidate chromosome equals an ancestor chromosome (gene_set == conf_set)
           → Remove the candidate (it's already in confirmed_ancestors)
        
        2. Nested: If candidate chromosome contains a complete ancestor chromosome 
           (conf_set is proper subset of gene_set)
           → Calculate difference as Residual with telomeres=True
           → Remove original chromosome
           → Create new Residual chromosome with range-based ID
        
        3. Neither: If neither Shared nor Nested (e.g., candidate is subset of ancestor, or no overlap at all)
           → Keep unchanged
        """
        
        total_removed = 0
        
        name_to_node = {}
        for node in self.tree.traverse():
            name_to_node[node.name] = node
        
        for node_name, nk in self.node_karyotypes.items():
            if nk.is_leaf:
                continue
            
            confirmed = self.confirmed_ancestors.get(node_name, {})
            
            confirmed_by_node = {}
            
            node_tree = name_to_node.get(node_name)
            if node_tree is not None:
                ancestor_list = list(node_tree.iter_ancestors())
                ancestor_list.reverse()
                for ancestor_node in ancestor_list:
                    anc_confirmed = self.confirmed_ancestors.get(ancestor_node.name, {})
                    if anc_confirmed:
                        confirmed_by_node[ancestor_node.name] = [
                            (cid, set(genes)) for cid, genes in anc_confirmed.items()
                        ]
            
            for cid, genes in confirmed.items():
                if node_name not in confirmed_by_node:
                    confirmed_by_node[node_name] = []
                confirmed_by_node[node_name].append((cid, set(genes)))
            
            confirmed_genes = set()
            for node_anc, chr_list in confirmed_by_node.items():
                for _, gene_set in chr_list:
                    confirmed_genes.update(gene_set)
            
            to_remove = []
            new_residuals = []

            for cid, genes in nk.chromosomes.items():
                gene_set = set(genes)
                matched = False
                
                for anc_node, chr_list in confirmed_by_node.items():
                    for conf_cid, conf_set in chr_list:
                        if gene_set == conf_set:
                            to_remove.append((cid, "equals_confirmed"))
                            matched = True
                            break
                        
                        if conf_set < gene_set:
                            remaining = gene_set - conf_set
                            remaining_ordered = [g for g in genes if g in remaining]
                            new_residuals.append((remaining_ordered, anc_node, cid))
                            to_remove.append((cid, "contains_ancestor->residual"))
                            matched = True
                            break
                    if matched:
                        break
            node_nk = self.node_karyotypes[node_name]
            current_node_confirmed_cids = set(confirmed.keys())  # CIDs of this node's confirmed ancestors
            
            for cid, reason in to_remove:
                if cid in current_node_confirmed_cids:
                    continue
                if cid in node_nk.chromosomes:
                    del node_nk.chromosomes[cid]
                if cid in node_nk.chr_attributes:
                    del node_nk.chr_attributes[cid]
                total_removed += 1
            
            if new_residuals:
                for genes, from_node, source_cid in new_residuals:
                    new_id = self._canonical_residual_id(source_cid, genes)
                    base_id = new_id
                    suffix = 2
                    while new_id in node_nk.chromosomes:
                        new_id = f"{base_id}_{suffix}"
                        suffix += 1
                    node_nk.chromosomes[new_id] = genes
                    node_nk.chr_attributes[new_id] = {
                        "provenance": f"Residual (from {from_node})",
                        "telomeres": True
                    }
        self.logs.append(f"--- Chromosome Cleanup Complete: removed {total_removed} chromosomes ---")

    def _index_confirmed_ancestors(self):
        allowed = {"Shared (Strict)", "Nested (Strict)"}
        confirmed = {}
        for name, nk in self.node_karyotypes.items():
            if nk.is_leaf:
                confirmed[name] = {}
                continue
            node_confirmed = {}
            for cid, genes in nk.chromosomes.items():
                attrs = nk.chr_attributes.get(cid, {})
                if attrs.get("telomeres", False) and attrs.get("provenance") in allowed:
                    node_confirmed[cid] = list(genes)
            by_set = {}
            for cid, genes in list(node_confirmed.items()):
                key = frozenset(genes)
                if key not in by_set:
                    by_set[key] = (cid, genes)
                    continue
                prev_cid, prev_genes = by_set[key]
                if len(prev_genes) != len(genes):
                    keep = prev_cid if len(prev_genes) < len(genes) else cid
                else:
                    keep = self._pick_preferred_equivalent_id(prev_cid, cid)
                drop = cid if keep == prev_cid else prev_cid
                keep_genes = prev_genes if keep == prev_cid else genes
                by_set[key] = (keep, keep_genes)
                if drop in node_confirmed:
                    del node_confirmed[drop]
                if drop in nk.chromosomes:
                    del nk.chromosomes[drop]
                if drop in nk.chr_attributes:
                    del nk.chr_attributes[drop]
            confirmed[name] = node_confirmed
        self.confirmed_ancestors = confirmed

        self._cleanup_redundant_chromosomes()

        self.logs.append("--- Confirmed Ancestors By Node ---")
        for node in self.tree.traverse("preorder"):
            if node.is_leaf():
                continue
            node_name = node.name
            node_label = "Root" if node.is_root() else node_name
            
            node_chrs = self.confirmed_ancestors.get(node_name, {})
            unique_chrs = {}
            seen_sets = set()
            for cid, genes in node_chrs.items():
                gset = frozenset(genes)
                if gset not in seen_sets:
                    seen_sets.add(gset)
                    unique_chrs[cid] = genes
            
            self.logs.append(f"[Confirmed@{node_label}] count={len(unique_chrs)}")
            nk = self.node_karyotypes.get(node_name)
            for cid, genes in unique_chrs.items():
                prov = ""
                telo = ""
                if nk is not None:
                    attrs = nk.chr_attributes.get(cid, {})
                    prov = attrs.get("provenance", "")
                    telo = attrs.get("telomeres", False)
                self.logs.append(f"  - {cid}\tgenes={len(genes)}\ttelomeres={telo}\tprovenance={prov}")
            self.logs.append("")

    def _apply_root_outgroup_decision_constraints(self):
        rows = getattr(self, "root_rct_outgroup_decisions", None) or []
        if not rows:
            return
        if not self.cfg.get("enable_root_outgroup_conflict_demote", True):
            return
        root = self.tree.get_tree_root()
        root_name = root.name
        nk_root = self.node_karyotypes.get(root_name)
        
        if not hasattr(self, "root_rct_survivor_chr_ids"):
            self.root_rct_survivor_chr_ids = set()
        if nk_root is None:
            return

        frozen_root_map = copy.deepcopy(self.frozen_confirmed_ancestors.get(root_name, {}) or {})
        frozen_root_sets = {frozenset(genes): cid for cid, genes in frozen_root_map.items()}

        def _weight(r):
            try:
                ls = int(r.get("LeftScore", 0))
                rs = int(r.get("RightScore", 0))
            except Exception:
                ls = 0
                rs = 0
            v = r.get("Verdict")
            if v == "left_is_ancestral":
                return ls - rs
            if v == "right_is_ancestral":
                return rs - ls
            return 0

        locked_proto = []
        locked_genes = set()
        root_rct_votes = {}

        def _pair_key(chr_ids):
            items = sorted([c for c in chr_ids if c])
            return "|".join(items)

        sorted_rows = sorted(rows, key=_weight, reverse=True)
        survivor_pair_keys = set()
        pair_source = {}
        for row in sorted_rows:
            verdict = row.get("Verdict")
            if verdict == "left_is_ancestral":
                anc_node = row.get("LeftNode")
                anc_chr_ids = [row.get("LeftChr1"), row.get("LeftChr2")]
                derived_chr_ids = [row.get("RightChr1"), row.get("RightChr2")]
            elif verdict == "right_is_ancestral":
                anc_node = row.get("RightNode")
                anc_chr_ids = [row.get("RightChr1"), row.get("RightChr2")]
                derived_chr_ids = [row.get("LeftChr1"), row.get("LeftChr2")]
            else:
                continue
            anc_pair_key = _pair_key(anc_chr_ids)
            derived_pair_key = _pair_key(derived_chr_ids)
            if anc_pair_key:
                survivor_pair_keys.add(anc_pair_key)
                pair_source[anc_pair_key] = (anc_node, list(anc_chr_ids))
            if derived_pair_key and derived_pair_key in survivor_pair_keys:
                survivor_pair_keys.remove(derived_pair_key)
        
        survivor_chr_ids = set()
        
        def _base_cid(x):
            s = str(x).split("#", 1)[0]
            if ":" not in s:
                return s
            left, right = s.split(":", 1)
            if "(" in left and ")" in left:
                return left
            if "(" in right and ")" in right:
                return right
            return left
        survivor_chr_bases = {_base_cid(cid) for cid in survivor_chr_ids}
        pruned_by_competition = 0
        overlap_skips = 0
        applied_rows = 0
        extracted_ops = 0
        pending_pairs = {}  # {pair_key: (anc_node, anc_chr_ids, proto_sets, row)} 保存因基因重叠被暂缓的胜者对

        for row in sorted_rows:
            verdict = row.get("Verdict")
            left_node = row.get("LeftNode")
            right_node = row.get("RightNode")
            if verdict == "left_is_ancestral":
                anc_node = left_node
                anc_chr_ids = [row.get("LeftChr1"), row.get("LeftChr2")]
                derived_node = right_node
                derived_chr_ids = [row.get("RightChr1"), row.get("RightChr2")]
            elif verdict == "right_is_ancestral":
                anc_node = right_node
                anc_chr_ids = [row.get("RightChr1"), row.get("RightChr2")]
                derived_node = left_node
                derived_chr_ids = [row.get("LeftChr1"), row.get("LeftChr2")]
            else:
                continue

            anc_pair_key = _pair_key(anc_chr_ids)
            derived_pair_key = _pair_key(derived_chr_ids)
            if anc_pair_key and anc_pair_key not in survivor_pair_keys:
                pruned_by_competition += 1
                continue
            winner_analysis_only = bool(row.get("WinnerAnalysisOnly", False))
            
            if winner_analysis_only:
                for cid in [c for c in anc_chr_ids if c]:
                    self.root_rct_survivor_chr_ids.add(cid)
                    if cid in self.confirmed_ancestors.get(root_name, {}):
                        pass
                    else:
                        genes = None
                        if cid in self.confirmed_ancestors.get(root_name, {}):
                            genes = self.confirmed_ancestors[root_name][cid]
                        elif cid in nk_root.chromosomes:
                            genes = nk_root.chromosomes[cid]
                        else:
                            for nk_name, nk in self.node_karyotypes.items():
                                if cid in nk.chromosomes:
                                    genes = nk.chromosomes[cid]
                                    break
                            if not genes:
                                for nk_name, node_chrs in self.confirmed_ancestors.items():
                                    if cid in node_chrs:
                                        genes = node_chrs[cid]
                                        break
                        if genes:
                            self.confirmed_ancestors.setdefault(root_name, {})[cid] = list(genes)
                continue
                
            for cid in [c for c in anc_chr_ids if c]:
                root_rct_votes[cid] = root_rct_votes.get(cid, 0) + 1
            for cid in [c for c in derived_chr_ids if c]:
                root_rct_votes[cid] = root_rct_votes.get(cid, 0) - 1
            anc_map = (self.confirmed_ancestors or {}).get(anc_node, {})
            nk_anc = self.node_karyotypes.get(anc_node)

            proto_sets = []
            for anc_cid in [c for c in anc_chr_ids if c]:
                target_genes = anc_map.get(anc_cid)
                if not target_genes and nk_anc is not None:
                    target_genes = nk_anc.chromosomes.get(anc_cid)
                if not target_genes and nk_root is not None:
                    target_genes = nk_root.chromosomes.get(anc_cid)
                if not target_genes:
                    for n_name, n_map in (self.confirmed_ancestors or {}).items():
                        if anc_cid in n_map:
                            target_genes = list(n_map.get(anc_cid) or [])
                            break
                if not target_genes:
                    for n_name, nk in (self.node_karyotypes or {}).items():
                        if anc_cid in nk.chromosomes:
                            target_genes = list(nk.chromosomes.get(anc_cid) or [])
                            break
                if not target_genes:
                    self.logs.append(f"[Root-Outgroup-Constraint]   {anc_cid}: genes not found anywhere")
                    continue
                target_set = set(target_genes)
                proto_sets.append((anc_cid, target_set))

            remaining_proto_sets = []
            for anc_cid, target_set in proto_sets:
                frozen_match_cid = frozen_root_sets.get(frozenset(target_set))
                if frozen_match_cid:
                    locked_proto.append((frozen_match_cid, target_set))
                    locked_genes |= target_set
                    survivor_chr_ids.add(frozen_match_cid)
                    self.confirmed_ancestors.setdefault(root_name, {})[frozen_match_cid] = list(frozen_root_map[frozen_match_cid])
                    self.logs.append(f"[Root-Outgroup-Constraint-FrozenKeep] {frozen_match_cid} already frozen; skip rewrite")
                    continue
                remaining_proto_sets.append((anc_cid, target_set))
            proto_sets = remaining_proto_sets
            if not proto_sets:
                continue

            conflict = False
            for _, pset in proto_sets:
                if pset & locked_genes:
                    conflict = True
                    break
            if conflict:
                pending_pairs[anc_pair_key] = (anc_node, anc_chr_ids, proto_sets, row)
                continue
            applied_rows += 1

            for anc_cid, target_set in proto_sets:
                target_genes = anc_map.get(anc_cid)
                if not target_genes and nk_anc is not None:
                    target_genes = nk_anc.chromosomes.get(anc_cid)
                if not target_genes:
                    target_genes = nk_root.chromosomes.get(anc_cid)
                if not target_genes:
                    continue
                candidates = []
                for rcid, rgenes in nk_root.chromosomes.items():
                    rset = set(rgenes)
                    if self._is_frozen_confirmed(root_name, cid=rcid, genes=rgenes) and rset != target_set:
                        continue
                    if target_set.issubset(rset):
                        candidates.append((len(rgenes), rcid))
                if not candidates:
                    continue
                candidates.sort()
                best_source = candidates[0][1]

                src_genes = nk_root.chromosomes.get(best_source, [])
                src_set = set(src_genes)
                if not target_set.issubset(src_set):
                    continue

                new_id = anc_cid
                if new_id in nk_root.chromosomes and set(nk_root.chromosomes.get(new_id, [])) != target_set:
                    new_id = f"{new_id}_root"

                residue_genes = [g for g in src_genes if g not in target_set]
                if best_source == new_id:
                    nk_root.chromosomes[new_id] = list(target_genes)
                    orig_attrs = nk_root.chr_attributes.get(best_source, {})
                    orig_prov = orig_attrs.get("provenance", "")
                    if orig_prov in {"Shared (Strict)", "Nested (Strict)"}:
                        nk_root.chr_attributes[new_id] = {"provenance": orig_prov, "telomeres": True}
                    else:
                        nk_root.chr_attributes[new_id] = {"provenance": "Reconciled (Strict-RCT)", "telomeres": True}
                    self.confirmed_ancestors.setdefault(root_name, {})[new_id] = list(target_genes)
                    if residue_genes:
                        residue_id = self._generate_range_id(f"{best_source}_residue", residue_genes)
                        nk_root.chromosomes[residue_id] = residue_genes
                        nk_root.chr_attributes[residue_id] = {"provenance": "Residual (RCT)", "telomeres": False}
                else:
                    nk_root.chromosomes[new_id] = list(target_genes)
                    src_attrs = nk_root.chr_attributes.get(best_source, {})
                    src_prov = src_attrs.get("provenance", "")
                    if src_prov in {"Shared (Strict)", "Nested (Strict)"}:
                        nk_root.chr_attributes[new_id] = {"provenance": src_prov, "telomeres": True}
                    else:
                        nk_root.chr_attributes[new_id] = {"provenance": "Reconciled (Strict-RCT)", "telomeres": True}
                    self.confirmed_ancestors.setdefault(root_name, {})[new_id] = list(target_genes)

                    if residue_genes:
                        residue_id = self._generate_range_id(f"{best_source}_residue", residue_genes)
                        nk_root.chromosomes[residue_id] = residue_genes
                        nk_root.chr_attributes[residue_id] = {"provenance": "Residual (RCT)", "telomeres": False}
                    if src_prov in {"Shared (Strict)", "Nested (Strict)"} or self._is_frozen_confirmed(root_name, cid=best_source, genes=src_genes):
                        extracted_ops += 1
                    else:
                        nk_root.chromosomes.pop(best_source, None)
                        nk_root.chr_attributes.pop(best_source, None)
                        if root_name in self.confirmed_ancestors:
                            self.confirmed_ancestors[root_name].pop(best_source, None)
                        extracted_ops += 1

                to_update = []
                for rcid, rgenes in list(nk_root.chromosomes.items()):
                    if rcid == new_id:
                        continue
                    rset = set(rgenes)
                    if target_set.issubset(rset):
                        to_update.append(rcid)
                for rcid in to_update:
                    rgenes = nk_root.chromosomes.get(rcid, [])
                    if self._is_frozen_confirmed(root_name, cid=rcid, genes=rgenes):
                        self.logs.append(f"[Root-Outgroup-Constraint-FrozenSkip] keep frozen source {rcid}")
                        continue
                    filtered = [g for g in rgenes if g not in target_set]
                    if not filtered:
                        attrs = nk_root.chr_attributes.get(rcid, {})
                        if attrs.get("provenance") in {"Shared (Strict)", "Nested (Strict)"}:
                            continue
                        nk_root.chromosomes.pop(rcid, None)
                        nk_root.chr_attributes.pop(rcid, None)
                        if root_name in self.confirmed_ancestors:
                            self.confirmed_ancestors[root_name].pop(rcid, None)
                        extracted_ops += 1
                        continue
                    nk_root.chromosomes[rcid] = filtered
                    if root_name in self.confirmed_ancestors and rcid in self.confirmed_ancestors.get(root_name, {}):
                        self.confirmed_ancestors[root_name][rcid] = filtered
                    extracted_ops += 1

            for cid, pset in proto_sets:
                locked_proto.append((cid, pset))
                locked_genes |= pset
                survivor_chr_ids.add(cid)  # 成功应用，添加到 survivor

            if derived_pair_key and derived_pair_key in pending_pairs:
                del pending_pairs[derived_pair_key]
            if anc_pair_key in pending_pairs:
                del pending_pairs[anc_pair_key]

            demoted = 0
            for rcid, rgenes in list(nk_root.chromosomes.items()):
                if rcid in [c for c, _ in locked_proto]:
                    continue
                rset = set(rgenes)
                overlaps = [p for p, pset in locked_proto if (rset & pset)]
                if len(overlaps) >= 2:
                    attrs = nk_root.chr_attributes.get(rcid, {})
                    if attrs.get("provenance") in {"Shared (Strict)", "Nested (Strict)"} or self._is_frozen_confirmed(root_name, cid=rcid, genes=rgenes):
                        continue
                    nk_root.chr_attributes[rcid] = {"provenance": "Derived (RCT-Conflict)", "telomeres": False}
                    if root_name in self.confirmed_ancestors:
                        self.confirmed_ancestors[root_name].pop(rcid, None)
                    demoted += 1

        pending_applied = 0
        pending_final_skips = 0
        for pair_key, (anc_node, anc_chr_ids, proto_sets, row) in list(pending_pairs.items()):
            if pair_key not in survivor_pair_keys:
                pending_final_skips += 1
                continue

            self.logs.append(f"[Root-Outgroup-Constraint] processing pending winner: {pair_key}")
            self.logs.append(f"[Root-Outgroup-Constraint]   proto_sets contains: {[c for c, _ in proto_sets]}")
            anc_map = (self.confirmed_ancestors or {}).get(anc_node, {})
            nk_anc = self.node_karyotypes.get(anc_node)
            pair_applied = False

            filtered_proto_sets = []
            for anc_cid, target_set in proto_sets:
                frozen_match_cid = frozen_root_sets.get(frozenset(target_set))
                if frozen_match_cid:
                    locked_proto.append((frozen_match_cid, target_set))
                    locked_genes |= target_set
                    survivor_chr_ids.add(frozen_match_cid)
                    self.confirmed_ancestors.setdefault(root_name, {})[frozen_match_cid] = list(frozen_root_map[frozen_match_cid])
                    self.logs.append(f"[Root-Outgroup-Constraint-FrozenKeep] pending exact frozen {frozen_match_cid}")
                    continue
                filtered_proto_sets.append((anc_cid, target_set))

            for anc_cid, target_set in filtered_proto_sets:
                if target_set & locked_genes:
                    self.logs.append(f"[Root-Outgroup-Constraint]   {anc_cid}: already applied (genes in locked_genes)")
                    continue

                applied_rows += 1
                pending_applied += 1
                pair_applied = True
                self.logs.append(f"[Root-Outgroup-Constraint]   applying {anc_cid}")

                target_genes = anc_map.get(anc_cid)
                if not target_genes and nk_anc is not None:
                    target_genes = nk_anc.chromosomes.get(anc_cid)
                if not target_genes and nk_root is not None:
                    target_genes = nk_root.chromosomes.get(anc_cid)
                if not target_genes:
                    continue
                candidates = []
                for rcid, rgenes in nk_root.chromosomes.items():
                    rset = set(rgenes)
                    if self._is_frozen_confirmed(root_name, cid=rcid, genes=rgenes) and rset != target_set:
                        continue
                    if target_set.issubset(rset):
                        candidates.append((len(rgenes), rcid))
                if not candidates:
                    continue
                candidates.sort()
                best_source = candidates[0][1]

                src_genes = nk_root.chromosomes.get(best_source, [])
                src_set = set(src_genes)
                if not target_set.issubset(src_set):
                    continue

                new_id = anc_cid
                if new_id in nk_root.chromosomes and set(nk_root.chromosomes.get(new_id, [])) != target_set:
                    new_id = f"{new_id}_root"

                residue_genes = [g for g in src_genes if g not in target_set]
                if best_source == new_id:
                    nk_root.chromosomes[new_id] = list(target_genes)
                    orig_attrs = nk_root.chr_attributes.get(best_source, {})
                    orig_prov = orig_attrs.get("provenance", "")
                    if orig_prov in {"Shared (Strict)", "Nested (Strict)"}:
                        nk_root.chr_attributes[new_id] = {"provenance": orig_prov, "telomeres": True}
                    else:
                        nk_root.chr_attributes[new_id] = {"provenance": "Reconciled (Strict-RCT-Pending)", "telomeres": True}
                    self.confirmed_ancestors.setdefault(root_name, {})[new_id] = list(target_genes)
                    if residue_genes:
                        residue_id = self._generate_range_id(f"{best_source}_residue", residue_genes)
                        nk_root.chromosomes[residue_id] = residue_genes
                        nk_root.chr_attributes[residue_id] = {"provenance": "Residual (RCT)", "telomeres": False}
                else:
                    nk_root.chromosomes[new_id] = list(target_genes)
                    src_attrs = nk_root.chr_attributes.get(best_source, {})
                    src_prov = src_attrs.get("provenance", "")
                    if src_prov in {"Shared (Strict)", "Nested (Strict)"}:
                        nk_root.chr_attributes[new_id] = {"provenance": src_prov, "telomeres": True}
                    else:
                        nk_root.chr_attributes[new_id] = {"provenance": "Reconciled (Strict-RCT-Pending)", "telomeres": True}
                    self.confirmed_ancestors.setdefault(root_name, {})[new_id] = list(target_genes)

                    if residue_genes:
                        residue_id = self._generate_range_id(f"{best_source}_residue", residue_genes)
                        nk_root.chromosomes[residue_id] = residue_genes
                        nk_root.chr_attributes[residue_id] = {"provenance": "Residual (RCT)", "telomeres": False}
                    if src_prov in {"Shared (Strict)", "Nested (Strict)"} or self._is_frozen_confirmed(root_name, cid=best_source, genes=src_genes):
                        extracted_ops += 1
                    else:
                        nk_root.chromosomes.pop(best_source, None)
                        nk_root.chr_attributes.pop(best_source, None)
                        if root_name in self.confirmed_ancestors:
                            self.confirmed_ancestors[root_name].pop(best_source, None)
                        extracted_ops += 1

                to_update = []
                for rcid, rgenes in list(nk_root.chromosomes.items()):
                    if rcid == new_id:
                        continue
                    rset = set(rgenes)
                    if target_set.issubset(rset):
                        to_update.append(rcid)
                for rcid in to_update:
                    rgenes = nk_root.chromosomes.get(rcid, [])
                    if self._is_frozen_confirmed(root_name, cid=rcid, genes=rgenes):
                        self.logs.append(f"[Root-Outgroup-Constraint-FrozenSkip] keep frozen source {rcid}")
                        continue
                    filtered = [g for g in rgenes if g not in target_set]
                    if not filtered:
                        attrs = nk_root.chr_attributes.get(rcid, {})
                        if attrs.get("provenance") in {"Shared (Strict)", "Nested (Strict)"}:
                            continue
                        nk_root.chromosomes.pop(rcid, None)
                        nk_root.chr_attributes.pop(rcid, None)
                        if root_name in self.confirmed_ancestors:
                            self.confirmed_ancestors[root_name].pop(rcid, None)
                        extracted_ops += 1
                        continue
                    nk_root.chromosomes[rcid] = filtered
                    if root_name in self.confirmed_ancestors and rcid in self.confirmed_ancestors.get(root_name, {}):
                        self.confirmed_ancestors[root_name][rcid] = filtered
                    extracted_ops += 1

                locked_proto.append((anc_cid, target_set))
                locked_genes |= target_set
                survivor_chr_ids.add(anc_cid)  # 成功应用，添加到 survivor

        if pending_pairs:
            self.logs.append(f"[Root-Outgroup-Constraint] pending_pairs processed: {len(pending_pairs)} total, {pending_applied} applied, {pending_final_skips} final skips")

        def _resolve_chr_genes(cid, preferred_node):
            p_map = (self.confirmed_ancestors or {}).get(preferred_node, {}) if preferred_node else {}
            if cid in p_map:
                return p_map.get(cid), f"{preferred_node}.confirmed"
            p_nk = self.node_karyotypes.get(preferred_node) if preferred_node else None
            if p_nk is not None and cid in p_nk.chromosomes:
                return p_nk.chromosomes.get(cid), f"{preferred_node}.karyotype"
            root_map = (self.confirmed_ancestors or {}).get(root_name, {})
            if cid in root_map:
                return root_map.get(cid), f"{root_name}.confirmed"
            if cid in nk_root.chromosomes:
                return nk_root.chromosomes.get(cid), f"{root_name}.karyotype"
            for n_name, n_map in (self.confirmed_ancestors or {}).items():
                if cid in n_map:
                    return n_map.get(cid), f"{n_name}.confirmed"
            for n_name, nk in (self.node_karyotypes or {}).items():
                if cid in nk.chromosomes:
                    return nk.chromosomes.get(cid), f"{n_name}.karyotype"
            return None, "missing"

        materialized_count = 0
        materialize_miss = 0
        for pair_key in sorted(survivor_pair_keys):
            anc_node, anc_chr_ids = pair_source.get(pair_key, (None, []))
            if not anc_node:
                continue
            for anc_cid in [c for c in anc_chr_ids if c]:
                if anc_cid in nk_root.chromosomes:
                    continue
                target_genes, source_tag = _resolve_chr_genes(anc_cid, anc_node)
                if not target_genes:
                    materialize_miss += 1
                    continue
                nk_root.chromosomes[anc_cid] = list(target_genes)
                nk_root.chr_attributes[anc_cid] = {"provenance": "Reconciled (Strict-RCT)", "telomeres": True}
                self.confirmed_ancestors.setdefault(root_name, {})[anc_cid] = list(target_genes)
                materialized_count += 1

        vote_demoted = 0
        for rcid in list(nk_root.chromosomes.keys()):
            vote = root_rct_votes.get(rcid, 0)
            if vote >= 0:
                continue
            if rcid in survivor_chr_ids or _base_cid(rcid) in survivor_chr_bases:
                continue
            attrs = nk_root.chr_attributes.get(rcid, {})
            if attrs.get("provenance") in {"Shared (Strict)", "Nested (Strict)"} or self._is_frozen_confirmed(root_name, cid=rcid, genes=nk_root.chromosomes.get(rcid, [])):
                continue
            nk_root.chr_attributes[rcid] = {"provenance": "Derived (RCT-Vote)", "telomeres": False}
            if root_name in self.confirmed_ancestors:
                self.confirmed_ancestors[root_name].pop(rcid, None)
            vote_demoted += 1
        survivor_chr_bases = {_base_cid(cid) for cid in survivor_chr_ids}
        self.root_rct_vote_map = dict(root_rct_votes)
        self.root_rct_survivor_chr_ids = set(survivor_chr_ids)
        self.root_rct_survivor_chr_bases = set(survivor_chr_bases)

    def _init_residual_pools(self):
        """Initialize residual pools from node_karyotypes.
        
        Residual chromosomes are chromosomes with telomeres=True that are NOT
        already confirmed in any ancestor node's confirmed_ancestors.
        
        Uses gene set comparison to handle cases where leaf node chromosome IDs
        differ from ancestor chromosome IDs.
        """
        self.resolved_chromosomes = {name: set() for name in self.node_karyotypes}
        self.residual_chromosomes = {}
        
        applicable_confirmed_cache = {}
        for name, nk in self.node_karyotypes.items():
            applicable_confirmed_cache[name] = self._get_applicable_confirmed_gene_sets(name, include_self=True)

        for name, nk in self.node_karyotypes.items():
            residual_pool = {}
            applicable_confirmed_gene_sets = applicable_confirmed_cache.get(name, set())
            for cid, genes in nk.chromosomes.items():
                attrs = nk.chr_attributes.get(cid, {})
                if not attrs.get("telomeres", False):
                    continue
                if frozenset(genes) in applicable_confirmed_gene_sets:
                    continue
                range_id = self._canonical_residual_id(cid, genes)
                residual_pool[range_id] = list(genes)
            self.residual_chromosomes[name] = residual_pool

    def _canonical_residual_id(self, cid, genes):
        if ":" in cid:
            return cid
        return self._generate_range_id(cid, genes)

    def _get_base_chr_id(self, cid):
        if ":" in cid:
            return cid.split(":")[0]
        return cid

    def _check_residual_telomeres(self, cid, genes, node_name, parent_confirmed, original_genes=None, parent_sets=None):
        real_cid = cid.split(":")[-1] if ":" in cid else cid
        leaf_name = cid.split(":")[0] if ":" in cid else node_name
        
        # 1. Check existing attributes in NodeKaryotype (e.g. from Leaf)
        for check_node in [leaf_name, node_name]:
            nk = self.node_karyotypes.get(check_node)
            if nk:
                attrs = nk.chr_attributes.get(real_cid, {})
                if attrs.get("telomeres", False):
                    full_genes = nk.chromosomes.get(real_cid, [])
                    if len(genes) == len(full_genes) and set(genes) == set(full_genes):
                        return True

        # 2. Check match with parent confirmed (use pre-computed sets if available)
        subset_set = set(genes)
        psets = parent_sets if parent_sets is not None else {pcid: set(pgenes) for pcid, pgenes in parent_confirmed.items()}
        for pcid, pset in psets.items():
            if pset == subset_set and len(parent_confirmed[pcid]) == len(genes):
                return True
        if original_genes:
            original_set = set(original_genes)
            remainder_set = original_set - subset_set
            if remainder_set:
                remainder_count = len(original_genes) - len(genes)
                for pcid, pset in psets.items():
                    if pset == remainder_set and len(parent_confirmed[pcid]) == remainder_count:
                        return True
        return False

    def _collect_branch_residuals(self, node, residual_local):
        """
        Collect residual chromosomes from all leaf nodes under the given node's branch.
        Returns a merged residual pool dictionary.
        """
        branch_pool = {}
        for leaf in node.iter_leaves():
            leaf_name = leaf.name
            if leaf_name in residual_local:
                for cid, genes in residual_local[leaf_name].items():
                    branch_pool[f"{leaf_name}:{cid}"] = list(genes)
        return branch_pool

    def _collect_sister_residuals(self, node, residual_local):
        """
        Collect residual chromosomes from all leaf nodes under the sister branch.
        Returns a merged residual pool dictionary.
        """
        if node.is_root():
            outgroup_k = getattr(self, "outgroup_karyotypes", None) or {}
            if not outgroup_k:
                return None
            sister_pool = {}
            for cid, genes in outgroup_k.items():
                sister_pool[f"Outgroup:{cid}"] = list(genes)
            return sister_pool

        sisters = node.get_sisters()
        if not sisters:
            return None

        sister_pool = {}
        for sister in sisters:
            for leaf in sister.iter_leaves():
                leaf_name = leaf.name
                if leaf_name in residual_local:
                    for cid, genes in residual_local[leaf_name].items():
                        sister_pool[f"{leaf_name}:{cid}"] = list(genes)
        return sister_pool if sister_pool else None

    def _has_nonidentical_overlap_with_sets(self, new_genes, gene_sets):
        """
        Return (True, overlap_size) if new_genes overlaps any set in gene_sets
        with a non-identical relationship. Exact duplicates are ignored here.
        """
        new_set = set(new_genes)
        for genes in gene_sets:
            other = set(genes)
            overlap = new_set & other
            if overlap and new_set != other:
                return True, len(overlap)
        return False, 0

    def _prune_overlapping_confirmed_for_node(self, node_name):
        """
        Enforce a hard invariant: confirmed ancestors at the same node must be pairwise disjoint.
        Keep earlier entries (dict insertion order) and demote later overlapping ones.
        """
        node_confirmed = self.confirmed_ancestors.get(node_name, {})
        if not node_confirmed:
            return
        nk = self.node_karyotypes.get(node_name)
        kept = {}
        occupied = set()
        dropped = []
        for cid, genes in node_confirmed.items():
            gset = set(genes)
            overlap = gset & occupied
            if overlap:
                dropped.append((cid, len(overlap), len(genes)))
                if nk is not None:
                    attrs = nk.chr_attributes.setdefault(cid, {})
                    old_prov = attrs.get('provenance', '')
                    attrs['provenance'] = (old_prov + ' | Overlap-Demoted').strip(' |')
                continue
            kept[cid] = genes
            occupied.update(gset)
        if dropped:
            self.confirmed_ancestors[node_name] = kept
            if nk is not None:
                for cid, _, _ in dropped:
                    nk.chr_attributes.setdefault(cid, {})['telomeres'] = False
            node_label = 'Root' if self.tree.get_tree_root().name == node_name else node_name
            for cid, ov, n in dropped:
                self.logs.append(f"[Confirmed-Overlap-Pruned] {node_label}: demoted {cid} genes={n} overlap={ov}")

    def _check_conflict_with_confirmed(self, node_name, new_genes):
        """
        Check if a new chromosome conflicts with confirmed ancestral chromosomes.
        For Root, any non-identical overlap is forbidden.
        For internal nodes, keep the older relaxed behavior for subset/superset relationships.
        Returns: (has_conflict, conflict_info)
        """
        confirmed = self.confirmed_ancestors.get(node_name, {})
        new_set = set(new_genes)
        root_name = self.tree.get_tree_root().name if self.tree is not None else None
        strict_root = (node_name == root_name)

        for cid, genes in confirmed.items():
            conf_set = set(genes)
            overlap = new_set & conf_set
            if not overlap:
                continue
            if new_set == conf_set:
                continue
            if strict_root:
                return True, f"overlap={len(overlap)} with {cid}"
            if not new_set < conf_set and not conf_set < new_set:
                return True, f"overlap={len(overlap)} with {cid}"

        return False, None

    def _discover_from_residuals(self, iteration=1):
        self.logs.append(f"--- Residual Discovery Pass (Nested/Shared on Residual Pools) [iter={iteration}] ---")
        residual_local = {name: copy.deepcopy(pool) for name, pool in self.residual_chromosomes.items()}
        added_any = False

        preferred_id_by_set = self._compute_preferred_ids()
        
        items_to_demote = []

        for node in self.tree.traverse("postorder"):
            if node.is_leaf():
                continue
            if len(node.children) < 2:
                continue
            node_label = node.name if node.name else "Root"
            n1 = node.children[0]
            n2 = node.children[1]
            n1_name = n1.name
            n2_name = n2.name
            
            pool1 = self._collect_branch_residuals(n1, residual_local)
            pool2 = self._collect_branch_residuals(n2, residual_local)
            
            current_pool = residual_local.get(node.name, {})
            if current_pool:
                child1_genes = set()
                child2_genes = set()
                for genes in pool1.values():
                    child1_genes.update(genes)
                for genes in pool2.values():
                    child2_genes.update(genes)
                
                for cid, genes in current_pool.items():
                    gene_set = set(genes)
                    overlap1 = len(gene_set & child1_genes)
                    overlap2 = len(gene_set & child2_genes)
                    
                    if overlap1 > overlap2:
                        pool1[cid] = genes
                    elif overlap2 > overlap1:
                        pool2[cid] = genes
            
            if not pool1 or not pool2:
                continue

            applicable_confirmed_sets = self._get_applicable_confirmed_gene_sets(node.name, include_self=True)

            stripped1 = {}
            stripped2 = {}
            for cid, genes in pool1.items():
                gene_set = frozenset(genes)
                if gene_set in applicable_confirmed_sets:
                    continue
                stripped1[cid] = list(genes)
            for cid, genes in pool2.items():
                gene_set = frozenset(genes)
                if gene_set in applicable_confirmed_sets:
                    continue
                stripped2[cid] = list(genes)

            parent_confirmed = self.confirmed_ancestors.get(node.name, {})
            min_block_size = int(self.cfg.get("min_block_size", 2))
            require_rct_support = bool(self.cfg.get("residual_rct_requires_support", True))

            parent_sets = {pcid: set(pgenes) for pcid, pgenes in parent_confirmed.items()}
            
            child1_confirmed = self.confirmed_ancestors.get(n1_name, {})
            child2_confirmed = self.confirmed_ancestors.get(n2_name, {})
            child1_sets = {pcid: set(pgenes) for pcid, pgenes in child1_confirmed.items()}
            child2_sets = {pcid: set(pgenes) for pcid, pgenes in child2_confirmed.items()}

            def has_cross_fusion(block_set, parent_set, chrom_genes):
                for i in range(len(chrom_genes) - 1):
                    a = chrom_genes[i]
                    b = chrom_genes[i + 1]
                    if (a in block_set and b in parent_set) or (a in parent_set and b in block_set):
                        return True, (a, b)
                return False, None

            def find_rct_support(block_set, child_pool, sibling_confirmed_sets, sibling_confirmed, sibling_name):
                """Find RCT support from parent ancestors AND sibling node's confirmed ancestors.
                
                Args:
                    block_set: The residual chromosome gene set to find support for
                    child_pool: The pool to search in (from the sibling node)
                    sibling_confirmed_sets: Gene sets of confirmed ancestors from sibling node
                    sibling_confirmed: Original confirmed ancestors dict from sibling node
                    sibling_name: Name of the sibling node (for logging)
                """
                best = None
                for ccid, cgenes in child_pool.items():
                    cset = set(cgenes)
                    block_hits = len(block_set & cset)
                    if block_hits < min_block_size:
                        continue
                    block_in_order = [g for g in cgenes if g in block_set]
                    if len(block_in_order) != block_hits:
                        continue
                    if not self._is_contiguous_block(block_in_order, cgenes, strict=True):
                        continue
                    
                    for pcid, pset in parent_sets.items():
                        if not pset.issubset(cset):
                            continue
                        parent_genes = parent_confirmed[pcid]
                        if not self._is_contiguous_block(parent_genes, cgenes, strict=True):
                            continue
                        parent_hits = len(pset)
                        ok, fp = has_cross_fusion(block_set, pset, cgenes)
                        if not ok:
                            continue
                        cand = (parent_hits, block_hits, ccid, pcid, fp, "parent")
                        if best is None or cand[:2] > best[:2]:
                            best = cand
                    
                    for scid, sset in sibling_confirmed_sets.items():
                        if not sset.issubset(cset):
                            continue
                        sibling_genes = sibling_confirmed[scid]
                        if not self._is_contiguous_block(sibling_genes, cgenes, strict=True):
                            continue
                        sibling_hits = len(sset)
                        ok, fp = has_cross_fusion(block_set, sset, cgenes)
                        if not ok:
                            continue
                        cand = (sibling_hits, block_hits, ccid, scid, fp, "sibling")
                        if best is None or cand[:2] > best[:2]:
                            best = cand
                
                return best

            existing_sets = {frozenset(v) for v in self.confirmed_ancestors.get(node.name, {}).values()}
            for ancestor in node.iter_ancestors():
                for v in self.confirmed_ancestors.get(ancestor.name, {}).values():
                    existing_sets.add(frozenset(v))

            enable_rct = self.cfg.get("enable_residual_rct_discovery", True) and (parent_confirmed or child1_confirmed or child2_confirmed)
            if enable_rct:
                def promote_intact_residue(source_name, source_cid, residue_genes, support, support_source="parent"):
                    nonlocal added_any
                    residue_set = frozenset(residue_genes)
                    if residue_set in existing_sets:
                        return
                    source_confirmed = self.confirmed_ancestors.get(source_name, {})
                    if source_cid in source_confirmed:
                        return
                    
                    parent_hits, block_hits, ccid, pcid, fp, support_type = support
                    support_source_label = f"{support_type}:{pcid}"
                    
                    require_outgroup = self.cfg.get("residual_rct_require_outgroup", True)
                    if require_outgroup:
                        outgroup_fp = self._get_outgroup_fusion_points(node)
                        if outgroup_fp is None:
                            self.logs.append(
                                f"[Residual-RCT] {node_label}: skip promote {source_name}:{source_cid} - no outgroup data"
                            )
                            return
                        normalized_fp = _normalize_pair(fp[0], fp[1])
                        if normalized_fp not in outgroup_fp:
                            self.logs.append(
                                f"[Residual-RCT] {node_label}: skip promote {source_name}:{source_cid} - "
                                f"fp={fp[0]}|{fp[1]} not in outgroup"
                            )
                            return
                    
                    new_id = self._canonical_residual_id(source_cid, residue_genes)
                    base_id = new_id
                    suffix = 2
                    while new_id in self.confirmed_ancestors.get(node.name, {}):
                        new_id = f"{base_id}_{suffix}"
                        suffix += 1
                    telomeres = self._check_residual_telomeres(source_cid, residue_genes, source_name, parent_confirmed, parent_sets=parent_sets)
                    self.logs.append(
                        f"[Residual-RCT] {node_label}: promote {new_id} from {source_name}:{source_cid} "
                        f"support={ccid} with {support_source_label} fp={fp[0]}|{fp[1]} hits={block_hits}/{parent_hits} telomeres={telomeres}"
                    )
                    nk = self.node_karyotypes.get(node.name)
                    if nk is not None:
                        nk.chromosomes[new_id] = list(residue_genes)
                        nk.chr_attributes[new_id] = {"provenance": "Residual (RCT)", "telomeres": telomeres}
                    if telomeres:
                        self.confirmed_ancestors.setdefault(node.name, {})[new_id] = list(residue_genes)
                    existing_sets.add(residue_set)
                    added_any = True

                if stripped1 and stripped2:
                    for cid, genes in stripped1.items():
                        if set(pool1.get(cid, [])) != set(genes):
                            continue
                        support = find_rct_support(set(genes), pool2, child2_sets, child2_confirmed, n2_name)
                        if support:
                            promote_intact_residue(n1_name, cid, genes, support)

                    for cid, genes in stripped2.items():
                        if set(pool2.get(cid, [])) != set(genes):
                            continue
                        support = find_rct_support(set(genes), pool1, child1_sets, child1_confirmed, n1_name)
                        if support:
                            promote_intact_residue(n2_name, cid, genes, support)

            if stripped1 and stripped2:
                index2 = {}
                for c2, genes2 in stripped2.items():
                    key = frozenset(genes2)
                    index2.setdefault(key, []).append(c2)

                used2 = set()
                for c1, genes1 in stripped1.items():
                    key = frozenset(genes1)
                    candidates = index2.get(key, [])
                    match = None
                    for c2 in candidates:
                        if c2 not in used2:
                            match = c2
                            break
                    if not match:
                        continue
                    used2.add(match)

                    rid1 = self._canonical_residual_id(c1, genes1)
                    rid2 = self._canonical_residual_id(match, stripped2[match])
                    shared_set = set(genes1)

                    support = None
                    if enable_rct:
                        support = find_rct_support(shared_set, pool1, child1_sets, child1_confirmed, n1_name) or find_rct_support(shared_set, pool2, child2_sets, child2_confirmed, n2_name)

                    prov = "RCT" if support is not None else "Shared"
                    telomeres = (set(pool1.get(c1, [])) == shared_set) or (set(pool2.get(match, [])) == shared_set)
                    self.logs.append(f"[Residual-{prov}] {node_label}: {n1_name}:{rid1} == {n2_name}:{rid2} (after-strip) telomeres={telomeres}")

                    new_set = frozenset(genes1)
                    if new_set not in existing_sets:
                        if self.cfg.get("enable_conflict_check", True):
                            has_conflict, conflict_info = self._check_conflict_with_confirmed(node.name, genes1)
                            if has_conflict:
                                self.logs.append(f"[Residual-{prov}-Conflict] {node_label}: rejected {n1_name}:{rid1} - {conflict_info}")
                                continue
                        
                        rid1_is_ancestor = ":" in rid1
                        rid2_is_ancestor = ":" in rid2
                        
                        new_id = None
                        if rid1_is_ancestor and rid2_is_ancestor:
                            if rid1 in self.confirmed_ancestors.get(node.name, {}):
                                new_id = rid1
                            elif rid2 in self.confirmed_ancestors.get(node.name, {}):
                                new_id = rid2
                            else:
                                new_id = rid1
                        elif rid1_is_ancestor:
                            new_id = rid1
                        elif rid2_is_ancestor:
                            new_id = rid2
                        else:
                            new_id = rid1
                        
                        nk = self.node_karyotypes.get(node.name)
                        if nk is not None:
                            nk.chromosomes[new_id] = list(genes1)
                            nk.chr_attributes[new_id] = {"provenance": f"Residual ({prov})", "telomeres": telomeres}
                        if telomeres:
                            self.confirmed_ancestors.setdefault(node.name, {})[new_id] = list(genes1)
                        added_any = True

            used1_nested = set()
            used2_nested = set()
            stripped1_sets = {c1: set(genes1) for c1, genes1 in stripped1.items()}
            stripped2_sets = {c2: set(genes2) for c2, genes2 in stripped2.items()}
            for c1, genes1 in stripped1.items():
                if c1 in used1_nested:
                    continue
                set1 = stripped1_sets[c1]
                for c2, genes2 in stripped2.items():
                    if c2 in used2_nested:
                        continue
                    set2 = stripped2_sets[c2]
                    if set1 < set2:
                        block = self._find_contiguous_block(genes2, set1)
                        if block:
                            rid1 = self._canonical_residual_id(c1, genes1)
                            rid2 = self._canonical_residual_id(c2, genes2)
                            telomeres = self._check_residual_telomeres(c1, genes1, n1_name, parent_confirmed, parent_sets=parent_sets)
                            if telomeres:
                                self.logs.append(f"[Residual-Nested] {node_label}: {n1_name}:{rid1} subset-of {n2_name}:{rid2} telomeres={telomeres}")
                            new_set = frozenset(genes1)
                            if new_set not in existing_sets:
                                if self.cfg.get("enable_conflict_check", True):
                                    has_conflict, conflict_info = self._check_conflict_with_confirmed(node.name, genes1)
                                    if has_conflict:
                                        self.logs.append(f"[Residual-Nested-Conflict] {node_label}: rejected {n1_name}:{rid1} - {conflict_info}")
                                        used1_nested.add(c1)
                                        used2_nested.add(c2)
                                        break
                                
                                new_id = rid1
                                nk = self.node_karyotypes.get(node.name)
                                if nk is not None:
                                    nk.chromosomes[new_id] = list(genes1)
                                    nk.chr_attributes[new_id] = {"provenance": "Residual (Nested)", "telomeres": telomeres}
                                if telomeres:
                                    self.confirmed_ancestors.setdefault(node.name, {})[new_id] = list(genes1)
                                existing_sets.add(new_set)
                                added_any = True
                            
                            telo2 = self._check_residual_telomeres(c2, genes2, n2_name, parent_confirmed, parent_sets=parent_sets)
                            
                            if self.cfg.get("enable_residual_nested_superset", True):
                                full_genes2 = pool2.get(c2, genes2)
                                superset_set = frozenset(full_genes2)
                                nk = self.node_karyotypes.get(node.name)
                                if nk is not None and not node.is_root():
                                    has_telo = self._check_residual_telomeres(c2, full_genes2, n2_name, parent_confirmed, parent_sets=parent_sets)
                                    if has_telo:
                                        confirmed_sets = [frozenset(g) for g in self.confirmed_ancestors.get(node.name, {}).values()]
                                        conflict_with = None
                                        for aset in confirmed_sets:
                                            if aset != superset_set and (aset & superset_set):
                                                conflict_with = aset
                                                break
                                        if conflict_with is not None:
                                            self.logs.append(f"[Residual-Nested-Superset-Conflict] {node_label}: disabled due to overlap={len(conflict_with & superset_set)} for {n2_name}:{c2}")
                                            self.cfg["enable_residual_nested_superset"] = False
                                            self._purge_residual_nested_supersets()
                                        else:
                                            preferred_id = c2
                                            nk.chromosomes[preferred_id] = list(full_genes2)
                                            nk.chr_attributes[preferred_id] = {"provenance": "Residual (Nested-Superset)", "telomeres": True}
                                            self.confirmed_ancestors.setdefault(node.name, {})[preferred_id] = list(full_genes2)
                                            existing_sets.add(superset_set)
                                            self.logs.append(f"[Residual-Nested-Superset] {node_label}: promoted {n2_name}:{preferred_id} (Telomeric Ancestor Form)")
                                            added_any = True

                            used1_nested.add(c1)
                            used2_nested.add(c2)
                            break
                    elif set2 < set1:
                        block = self._find_contiguous_block(genes1, set2)
                        if block:
                            rid1 = self._canonical_residual_id(c1, genes1)
                            rid2 = self._canonical_residual_id(c2, genes2)
                            telomeres = self._check_residual_telomeres(c2, genes2, n2_name, parent_confirmed, parent_sets=parent_sets)
                            if telomeres:
                                self.logs.append(f"[Residual-Nested] {node_label}: {n2_name}:{rid2} subset-of {n1_name}:{rid1} telomeres={telomeres}")
                            new_set = frozenset(genes2)
                            if new_set not in existing_sets:
                                if self.cfg.get("enable_conflict_check", True):
                                    has_conflict, conflict_info = self._check_conflict_with_confirmed(node.name, genes2)
                                    if has_conflict:
                                        self.logs.append(f"[Residual-Nested-Conflict] {node_label}: rejected {n2_name}:{rid2} - {conflict_info}")
                                        used1_nested.add(c1)
                                        used2_nested.add(c2)
                                        break
                                
                                new_id = rid2
                                nk = self.node_karyotypes.get(node.name)
                                if nk is not None:
                                    nk.chromosomes[new_id] = list(genes2)
                                    nk.chr_attributes[new_id] = {"provenance": "Residual (Nested)", "telomeres": telomeres}
                                if telomeres:
                                    self.confirmed_ancestors.setdefault(node.name, {})[new_id] = list(genes2)
                                existing_sets.add(new_set)
                                added_any = True
                            
                            if self.cfg.get("enable_residual_nested_superset", True):
                                full_genes1 = pool1.get(c1, genes1)
                                superset_set = frozenset(full_genes1)
                                nk = self.node_karyotypes.get(node.name)
                                if nk is not None and not node.is_root():
                                    has_telo = self._check_residual_telomeres(c1, full_genes1, n1_name, parent_confirmed, parent_sets=parent_sets)
                                    if has_telo:
                                        confirmed_sets = [frozenset(g) for g in self.confirmed_ancestors.get(node.name, {}).values()]
                                        conflict_with = None
                                        for aset in confirmed_sets:
                                            if aset != superset_set and (aset & superset_set):
                                                conflict_with = aset
                                                break
                                        if conflict_with is not None:
                                            self.logs.append(f"[Residual-Nested-Superset-Conflict] {node_label}: disabled due to overlap={len(conflict_with & superset_set)} for {n1_name}:{c1}")
                                            self.cfg["enable_residual_nested_superset"] = False
                                            self._purge_residual_nested_supersets()
                                        else:
                                            preferred_id = c1
                                            nk.chromosomes[preferred_id] = list(full_genes1)
                                            nk.chr_attributes[preferred_id] = {"provenance": "Residual (Nested-Superset)", "telomeres": True}
                                            self.confirmed_ancestors.setdefault(node.name, {})[preferred_id] = list(full_genes1)
                                            existing_sets.add(superset_set)
                                            self.logs.append(f"[Residual-Nested-Superset] {node_label}: promoted {n1_name}:{preferred_id} (Telomeric Ancestor Form)")
                                            added_any = True
                                    else:
                                        self.logs.append(f"[Residual-Nested-Superset-Rejected] {node_label}: {n1_name}:{rid1} has_telo=False")

                            used1_nested.add(c1)
                            used2_nested.add(c2)
                            break

        for node_name, pcid, node_label in items_to_demote:
            if self._is_frozen_confirmed(node_name, cid=pcid):
                self.logs.append(f"[Residual-RCT-Confirmed-Demoted-SkipFrozen] {node_label}: keep frozen {pcid}")
                continue
            if node_name in self.confirmed_ancestors and pcid in self.confirmed_ancestors[node_name]:
                del self.confirmed_ancestors[node_name][pcid]
                self.logs.append(f"[Residual-RCT-Confirmed-Demoted] {node_label}: demoted {pcid} from confirmed ancestors (endpoint check failed, kept in karyotypes for future iterations)")

        return added_any

    def _find_contiguous_block(self, genes, subset):
        subset_set = set(subset)
        block_start = None
        for i, g in enumerate(genes):
            if g in subset_set:
                if block_start is None:
                    block_start = i
            else:
                if block_start is not None:
                    block = genes[block_start:i]
                    if set(block) == subset_set:
                        return block
                    block_start = None
        if block_start is not None:
            block = genes[block_start:]
            if set(block) == subset_set:
                return block
        return None

    def _merge_redundant_chromosomes(self, chromosomes, attributes):
        pool = copy.deepcopy(chromosomes)
        pool_attrs = copy.deepcopy(attributes)
        
        while True:
            merged = False
            ids = sorted(pool.keys(), key=lambda k: len(pool[k]), reverse=True)
            
            for i in range(len(ids)):
                if merged: break
                id1 = ids[i]
                genes1 = set(pool[id1])
                
                for j in range(i + 1, len(ids)):
                    id2 = ids[j]
                    genes2 = set(pool[id2])
                    
                    overlap = len(genes1 & genes2)
                    if overlap == 0: continue
                    
                    jaccard = overlap / len(genes1 | genes2)
                    t1 = pool_attrs.get(id1, {}).get("telomeres", False)
                    t2 = pool_attrs.get(id2, {}).get("telomeres", False)
                    if t1 and t2 and jaccard < 0.999999:
                        continue
                    
                    if jaccard > 0.5:
                        pool[id1] = list(genes1 | genes2)
                        
                        attr1 = pool_attrs.get(id1, {})
                        attr2 = pool_attrs.get(id2, {})
                        pool_attrs[id1]["provenance"] = "Merged"
                        pool_attrs[id1]["telomeres"] = attr1.get("telomeres", False) or attr2.get("telomeres", False)
                        
                        del pool[id2]
                        if id2 in pool_attrs: del pool_attrs[id2]
                        
                        merged = True
                        break
            
            if not merged:
                break
                
        return pool, pool_attrs

    def _infer_internal_node(self, node):
        children = node.children
        nk = NodeKaryotype(node.name)
        node_label = node.name if node.name else "Root"
        
        if len(children) < 2:
            child_karyo = self.node_karyotypes[children[0].name]
            nk.chromosomes = copy.deepcopy(child_karyo.chromosomes)
            nk.chr_attributes = copy.deepcopy(child_karyo.chr_attributes)
            nk.add_log(f"Inherited directly from {children[0].name}")
        else:
            ordered_children = list(children)
            ordered_children.sort(key=lambda n: str(n.name))

            acc_karyo = self.node_karyotypes[ordered_children[0].name]
            acc_label = ordered_children[0].name
            acc_node = acc_karyo
            for child in ordered_children[1:]:
                right_karyo = self.node_karyotypes[child.name]
                nk.inference_log.append(f"--- Merge @Node {node_label}: left={acc_label} right={child.name} ---")
                inferred_chrs, attributes, logs = self._compare_and_merge(acc_node, right_karyo, node_label)
                nk.inference_log.extend(logs)
                tmp = NodeKaryotype(node.name)
                tmp.chromosomes = inferred_chrs
                tmp.chr_attributes = copy.deepcopy(acc_node.chr_attributes)
                tmp.chr_attributes.update(attributes)
                acc_node = tmp
                acc_label = node_label

            nk.chromosomes = acc_node.chromosomes
            nk.chr_attributes = acc_node.chr_attributes
            
        self.node_karyotypes[node.name] = nk
        self.logs.append(f"--- Node {node_label} Inference ---")
        self.logs.extend(nk.inference_log)
        self.logs.append("")

    def _compute_chr_telomeres(self, cid, genes):
        """Return True only if genes cover both ends of the source chromosome."""
        if ':' not in cid:
            return True
        base = cid.split(':', 1)[0]
        sp_prefix = base.split('(')[0]
        ref_genes = self.species_karyotypes.get(sp_prefix, {}).get(base)
        if not ref_genes:
            return True
        gene_set = frozenset(genes)
        return (ref_genes[0] in gene_set) and (ref_genes[-1] in gene_set)

    def _generate_range_id(self, base_id, genes):
        if not genes: return base_id
        clean_base = base_id.split(':', 1)[0].split('#', 1)[0]
        match = re.match(r"([^(]+)\(([^)]+)\)", clean_base)
        if not match:
            return f"{clean_base}:UnknownRange"
        sp_name = match.group(1)
        
        indices = []
        gene_set = set(genes)
        ref_genes = self.species_karyotypes.get(sp_name, {}).get(clean_base)
        if ref_genes:
            pos_in_chr = {}
            for i, g in enumerate(ref_genes, start=1):
                if g not in pos_in_chr:
                    pos_in_chr[g] = i
            hits = 0
            for g in gene_set:
                p = pos_in_chr.get(g)
                if p is None:
                    continue
                hits += 1
                indices.append(p)
            if hits < max(1, int(len(gene_set) * 0.8)):
                indices = []

        if not indices and sp_name in self.gene_position_map:
            sp_map = self.gene_position_map[sp_name]
            for g in genes:
                if g in sp_map:
                    indices.append(sp_map[g][1])
        
        if not indices: return f"{clean_base}:UnknownRange"
        indices = sorted(set(indices))
        
        if ref_genes:
            ref_set = set(ref_genes)
            if gene_set == ref_set:
                return clean_base  # Whole chromosome, no range needed
        
        blocks = []
        if indices:
            current_block = [indices[0]]
            for i in range(1, len(indices)):
                if indices[i] == indices[i-1] + 1:
                    current_block.append(indices[i])
                else:
                    blocks.append(current_block)
                    current_block = [indices[i]]
            blocks.append(current_block)
            
        range_strs = [f"{b[0]}-{b[-1]}" if b[0] != b[-1] else f"{b[0]}" for b in blocks]
        return f"{clean_base}:{','.join(range_strs)}"

    def _parse_sp_chr_from_cid(self, cid):
        """Parse species prefix and chr id from chromosome id like 'Sp1(1)'."""
        base = cid.split("#", 1)[0].split(":", 1)[0]
        m = re.match(r"([^(]+)\(([^)]+)\)", base)
        if not m:
            return None, None
        return m.group(1), m.group(2)

    def _sp_prefix_from_cid(self, cid):
        """Extract species prefix from chromosome id."""
        sp, _ = self._parse_sp_chr_from_cid(cid)
        if sp:
            return sp
        return cid.split("#", 1)[0].split(":", 1)[0]

    def _compute_preferred_ids(self):
        """Compute preferred chromosome IDs for each gene set across all confirmed ancestors."""
        id_counts_by_set = {}
        for _, node_chrs in self.confirmed_ancestors.items():
            for cid, genes in node_chrs.items():
                key = frozenset(genes)
                if key not in id_counts_by_set:
                    id_counts_by_set[key] = {}
                id_counts_by_set[key][cid] = id_counts_by_set[key].get(cid, 0) + 1

        preferred_id_by_set = {}
        for key, counts in id_counts_by_set.items():
            best = None
            for cid, cnt in counts.items():
                is_pending = "(pending)" in cid
                candidate = (0 if is_pending else 1, cnt, -len(cid), cid)
                if best is None or candidate > best:
                    best = candidate
            preferred_id_by_set[key] = best[3] if best is not None else next(iter(counts.keys()))
        return preferred_id_by_set

    def _cid_base(self, cid):
        return cid.split("#", 1)[0].split(":", 1)[0]

    def _ref_chr_len(self, cid):
        base = self._cid_base(cid)
        m = re.match(r"([^(]+)\(([^)]+)\)", base)
        if not m:
            return None
        sp = m.group(1)
        ref = self.species_karyotypes.get(sp, {}).get(base)
        if not ref:
            return None
        return len(ref)

    def _is_full_reference_chromosome(self, cid, genes_set):
        base = self._cid_base(cid)
        m = re.match(r"([^(]+)\(([^)]+)\)", base)
        if not m:
            return False
        sp = m.group(1)
        ref = self.species_karyotypes.get(sp, {}).get(base)
        if not ref:
            return False
        return len(ref) == len(genes_set) and set(ref) == genes_set

    def _pick_preferred_equivalent_id(self, cid1, cid2):
        return cid1 if len(cid1) <= len(cid2) else cid2

    def _get_outgroup_fusion_points(self, node):
        """
        Get fusion points from the node's sister branch (Outgroup).
        Unified logic for all nodes: use sister branch chromosomes as Outgroup.
        - Root node: sister branch = external outgroup_karyotypes
        - Non-Root node: sister branch = nodes on the other side of the tree and their descendants
        """
        if node.is_root():
            outgroup_k = getattr(self, "outgroup_karyotypes", None) or {}
            if not outgroup_k:
                return None
            return extract_fusion_points(outgroup_k)
        
        sisters = node.get_sisters()
        if not sisters:
            return None
        
        sister_karyotypes = {}
        for sister in sisters:
            sister_name = sister.name
            nk = self.node_karyotypes.get(sister_name)
            if nk:
                for cid, genes in nk.chromosomes.items():
                    sister_karyotypes[f"{sister_name}:{cid}"] = genes
            ca = self.confirmed_ancestors.get(sister_name, {})
            for cid, genes in ca.items():
                sister_karyotypes[f"{sister_name}:{cid}"] = genes
            for desc in sister.iter_descendants():
                desc_name = desc.name
                nk = self.node_karyotypes.get(desc_name)
                if nk:
                    for cid, genes in nk.chromosomes.items():
                        sister_karyotypes[f"{desc_name}:{cid}"] = genes
                ca = self.confirmed_ancestors.get(desc_name, {})
                for cid, genes in ca.items():
                    sister_karyotypes[f"{desc_name}:{cid}"] = genes
        
        if not sister_karyotypes:
            return None
        return extract_fusion_points(sister_karyotypes)

    def _is_contiguous_block(self, small_genes, large_genes, strict=False):
        """Check if small_genes forms a contiguous block within large_genes.
        
        Uses optimized set operations and list comprehension for performance.
        """
        if not small_genes:
            return False
        
        large_map = {g: i for i, g in enumerate(large_genes)}
        
        indices = [large_map[g] for g in small_genes if g in large_map]
        
        if not indices:
            return False
        
        indices.sort()
        n = len(indices)
        
        if strict:
            for i in range(n - 1):
                if indices[i + 1] != indices[i] + 1:
                    return False
            return True
        else:
            if n / len(small_genes) < 0.95:
                return False
            gaps = sum(1 for i in range(n - 1) if indices[i + 1] - indices[i] > 3)
            return gaps <= 2

    def _compare_and_merge(self, k1, k2, node_name):
        inferred = {}
        attributes = {}
        logs = []
        
        pool_k1 = {cid: list(genes) for cid, genes in k1.chromosomes.items()}
        pool_k2 = {cid: list(genes) for cid, genes in k2.chromosomes.items()}
        pool_attrs_k1 = {} 
        pool_attrs_k2 = {}
        
        def get_telo_k1(cid):
            return pool_attrs_k1.get(cid, k1.chr_attributes.get(cid, {}).get("telomeres", False))
            
        def get_telo_k2(cid):
            return pool_attrs_k2.get(cid, k2.chr_attributes.get(cid, {}).get("telomeres", False))

        def peel_residue(target_pool, target_attrs, target_id, residue_genes, prefix):
            if len(residue_genes) >= self.cfg["min_block_size"]:
                orig = target_pool[target_id]
                residue_ordered = [g for g in orig if g in residue_genes]
                if not residue_ordered: residue_ordered = list(residue_genes)
                
                new_id = self._generate_range_id(target_id, residue_ordered)
                target_pool[new_id] = residue_ordered
                target_attrs[new_id] = False # Fragment
                idxs = [i for i, g in enumerate(orig) if g in residue_genes]
                end_attached = False if not idxs else (min(idxs) == 0 or max(idxs) == len(orig) - 1)
                return (new_id, len(residue_ordered), end_attached)
            return None

        def conflicts_with_current_inferred(genes):
            for existing_id, existing_genes in inferred.items():
                overlap = set(genes) & set(existing_genes)
                if overlap and set(genes) != set(existing_genes):
                    return True, existing_id, len(overlap)
            return False, None, 0

        child_node_names = {k1.name, k2.name}
        ancestor_node_names = set()
        _cur = self._get_node_by_name(node_name)
        if _cur is not None:
            for _anc in _cur.iter_ancestors():
                ancestor_node_names.add(_anc.name)

        pass1_changed = True
        while pass1_changed:
            pass1_changed = False

            # --- 1.1 Strict Shared (Jaccard=1.0) ---
            candidates = []
            skip_k1 = set()
            skip_k2 = set()
            for id1 in list(pool_k1):
                if not get_telo_k1(id1): continue
                genes1 = set(pool_k1[id1])
                for id2 in list(pool_k2):
                    if not get_telo_k2(id2): continue
                    genes2 = set(pool_k2[id2])

                    if len(genes1) != len(genes2): continue

                    if genes1 == genes2:
                        candidates.append((1.0, 0, id1, id2))
            
            if candidates:
                for _, _, id1, id2 in candidates:
                    if id1 in pool_k1 and id2 in pool_k2:
                        out_id = self._pick_preferred_equivalent_id(id1, id2)
                        out_genes = list(pool_k1[id1]) if out_id == id1 else list(pool_k2[id2])
                        has_local_conflict, conflict_id, overlap_n = conflicts_with_current_inferred(out_genes)
                        if has_local_conflict:
                            logs.append(f"[Shared-Strict-LocalConflict] {id1} == {id2} rejected due to overlap={overlap_n} with inferred:{conflict_id}")
                            continue
                        inferred[out_id] = out_genes
                        attributes[out_id] = {"provenance": "Shared (Strict)", "telomeres": True}
                        logs.append(f"[Shared-Strict] {id1} == {id2} -> {out_id}")
                        del pool_k1[id1]
                        del pool_k2[id2]
                        # If a range-based fragment was matched, trim its base chromosome
                        # in the same pool to prevent it from being reused as a nested
                        # container with already-claimed genes.
                        matched_set = frozenset(out_genes)
                        if ":" in id1:
                            base1 = id1.split(":")[0]
                            if base1 in pool_k1:
                                trimmed1 = [g for g in pool_k1[base1] if g not in matched_set]
                                if trimmed1:
                                    pool_k1[base1] = trimmed1
                                else:
                                    del pool_k1[base1]
                                    pool_attrs_k1.pop(base1, None)
                        if ":" in id2:
                            base2 = id2.split(":")[0]
                            if base2 in pool_k2:
                                trimmed2 = [g for g in pool_k2[base2] if g not in matched_set]
                                if trimmed2:
                                    pool_k2[base2] = trimmed2
                                else:
                                    del pool_k2[base2]
                                    pool_attrs_k2.pop(base2, None)
                        pass1_changed = True
            
            if pass1_changed: continue 

            # --- 1.2 Strict Nested (100% Coverage, Contiguous, No Gaps) ---
            nested_changed = False
            for id1 in sorted(pool_k1.keys(), key=lambda k: len(pool_k1[k]), reverse=True):
                if id1 not in pool_k1: continue
                if not get_telo_k1(id1): continue
                genes1 = pool_k1[id1]

                best_target = None
                for id2 in pool_k2:
                    if not get_telo_k2(id2): continue
                    if len(pool_k2[id2]) >= len(genes1):
                        if set(genes1).issubset(set(pool_k2[id2])):
                            if self._is_contiguous_block(genes1, pool_k2[id2], strict=True):
                                best_target = id2
                                break
                
                if best_target:
                    id2 = best_target
                    genes2 = pool_k2[id2]
                    set1 = set(genes1)
                    set2 = set(genes2)
                    out_id = id1
                    out_genes = genes1
                    peeled_result = None
                    is_shared = False
                    if set1 == set2:
                        out_id = self._pick_preferred_equivalent_id(id1, id2)
                        out_genes = list(genes1) if out_id == id1 else list(genes2)
                        is_shared = True
                    else:
                        residue = set2 - set1
                        peeled_result = peel_residue(pool_k2, pool_attrs_k2, id2, residue, "NestedResidue")
                        if peeled_result is None:
                            logs.append(f"[Nested-Strict-PeelFail] {id1} in {id2}: peel failed, skipping")
                            continue
                    has_local_conflict, conflict_id, overlap_n = conflicts_with_current_inferred(out_genes)
                    if has_local_conflict:
                        logs.append(f"[Nested-Strict-LocalConflict] {id1} in {id2} rejected due to overlap={overlap_n} with inferred:{conflict_id}")
                        continue
                    inferred[out_id] = out_genes
                    if is_shared:
                        attributes[out_id] = {"provenance": "Shared (Strict)", "telomeres": True}
                        logs.append(f"[Shared-Strict] {id1} == {id2} -> {out_id}")
                    else:
                        attributes[out_id] = {"provenance": "Nested (Strict)", "telomeres": True}
                        peeled_id = peeled_result[0] if peeled_result else None
                        logs.append(f"[Nested-Strict] {id1} in {id2} -> {out_id} + Residue {peeled_id}")
                        if peeled_result:
                            logs.append(f"[Residue] {id2} -> {peeled_result[0]} size={peeled_result[1]} telomeres={peeled_result[2]}")
                    del pool_k1[id1]
                    del pool_k2[id2]
                    nested_changed = True
                    pass1_changed = True
                    break 
            
            if nested_changed: continue
            
            for id2 in sorted(pool_k2.keys(), key=lambda k: len(pool_k2[k]), reverse=True):
                if id2 not in pool_k2: continue
                if not get_telo_k2(id2): continue
                genes2 = pool_k2[id2]

                best_target = None
                for id1 in pool_k1:
                    if not get_telo_k1(id1): continue
                    if len(pool_k1[id1]) >= len(genes2):
                        if set(genes2).issubset(set(pool_k1[id1])):
                            if self._is_contiguous_block(genes2, pool_k1[id1], strict=True):
                                best_target = id1
                                break
                            
                if best_target:
                    id1 = best_target
                    genes1 = pool_k1[id1]
                    set2 = set(genes2)
                    set1 = set(genes1)
                    out_id = id2
                    out_genes = genes2
                    peeled_result = None
                    is_shared = False
                    if set1 == set2:
                        out_id = self._pick_preferred_equivalent_id(id1, id2)
                        out_genes = list(genes2) if out_id == id2 else list(genes1)
                        is_shared = True
                    else:
                        residue = set1 - set2
                        peeled_result = peel_residue(pool_k1, pool_attrs_k1, id1, residue, "NestedResidue")
                        if peeled_result is None:
                            logs.append(f"[Nested-Strict-PeelFail] {id2} in {id1}: peel failed, skipping")
                            continue
                    inferred[out_id] = out_genes
                    if is_shared:
                        attributes[out_id] = {"provenance": "Shared (Strict)", "telomeres": True}
                        logs.append(f"[Shared-Strict] {id2} == {id1} -> {out_id}")
                    else:
                        attributes[out_id] = {"provenance": "Nested (Strict)", "telomeres": True}
                        peeled_id = peeled_result[0] if peeled_result else None
                        logs.append(f"[Nested-Strict] {id2} in {id1} -> {out_id} + Residue {peeled_id}")
                        if peeled_result:
                            logs.append(f"[Residue] {id1} -> {peeled_result[0]} size={peeled_result[1]} telomeres={peeled_result[2]}")
                    del pool_k2[id2]
                    del pool_k1[id1]
                    nested_changed = True
                    pass1_changed = True
                    break
            
            if nested_changed: continue

        # ---------------- UNMATCHED ----------------
        for id1, genes in pool_k1.items():
            if id1 not in inferred:
                inferred[id1] = genes
                telo1 = get_telo_k1(id1)
                attributes[id1] = {"provenance": "Residue" if telo1 else "Unmatched", "telomeres": telo1}

        for id2, genes in pool_k2.items():
            if id2 not in inferred:
                inferred[id2] = genes
                telo2 = get_telo_k2(id2)
                attributes[id2] = {"provenance": "Residue" if telo2 else "Unmatched", "telomeres": telo2}
        return inferred, attributes, logs

# -----------------------------
# 4. VISUALIZATION
# -----------------------------
class TreeVisualizer:
    def __init__(self, tree):
        self.tree = tree
        
    def visualize_with_karyotypes(self, node_karyotypes, out_path):
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.scale = 80
        
        for node in self.tree.traverse():
            if not node.is_leaf():
                node.add_face(TextFace(node.name, fsize=12, fgcolor="red", bold=True), column=0, position="branch-top")

        self.tree.render(out_path, tree_style=ts)
        print(f"Tree visualization saved to {out_path}")

# -----------------------------
# 5. MAIN
# -----------------------------
def run(config=None):
    cfg = config if config is not None else CONFIG
    os.makedirs(cfg["output_dir"], exist_ok=True)
    loader = DataLoader(cfg)
    tree, sp_karyos, gene_pos_map = loader.load_data()
    
    fq_sp_karyos = {sp: {f"{sp}({cid})": genes for cid, genes in chrs.items()} for sp, chrs in sp_karyos.items()}
    
    reconstructor = AncestorReconstructor(tree, fq_sp_karyos, gene_pos_map, cfg=cfg)
    reconstructor.outgroup_name = cfg.get("outgroup_name", "Outgroup")
    reconstructor.outgroup_karyotypes = dict(getattr(loader, "outgroup_karyotypes", {}) or {})
    ancestors = reconstructor.run()

    decision_rows = getattr(reconstructor, "root_rct_outgroup_decisions", None)
    if decision_rows:
        out_path = os.path.join(cfg["output_dir"], "root_rct_outgroup_decisions.tsv")
        keys = [
            "ComparisonID",
            "LeftNode",
            "LeftChr1",
            "LeftChr2",
            "RightNode",
            "RightChr1",
            "RightChr2",
            "LeftScore",
            "RightScore",
            "Verdict",
            "LeftOrigPairs",
            "LeftTransPairs",
            "LeftOutgroupHit",
            "RightOrigPairs",
            "RightTransPairs",
            "RightOutgroupHit",
        ]
        with open(out_path, "w", encoding="utf-8") as f:
            f.write("\t".join(keys) + "\n")
            for row in decision_rows:
                f.write("\t".join(str(row.get(k, "")) for k in keys) + "\n")
    
    # Note: Detailed gene repertoire is now in gene_repertoire.txt
    with open(os.path.join(cfg["output_dir"], "reconstructed_ancestors.txt"), "w") as f:
        for anc in ancestors:
            f.write(f"ID: {anc['id']}\nGenes: {' '.join(anc['genes'])}\n{'-'*20}\n")
            
    with open(os.path.join(cfg["output_dir"], "inference_log.txt"), "w") as f:
        for line in reconstructor.logs:
            f.write(line + "\n")
            
    if cfg.get("enable_tree_viz", True):
        viz = TreeVisualizer(tree)
        tree_viz_dir = cfg.get("tree_viz_output_dir", cfg.get("output_dir", "output_reconstruction"))
        os.makedirs(tree_viz_dir, exist_ok=True)
        viz.visualize_with_karyotypes(reconstructor.node_karyotypes, os.path.join(tree_viz_dir, "tree_with_ancestors.png"))

    if cfg.get("enable_root_dotplot", True):
        root_tsv_path = os.path.join(cfg["output_dir"], "ancestor_gene_sets_by_node.tsv")
        true_root_path = cfg.get("input_true_root_karyotype", "")
        out_path = os.path.join(cfg["output_dir"], "DotPlot_ReconstructedRoot_vs_TrueRoot.png")
        generate_root_vs_trueroot_dotplot(root_tsv_path, true_root_path, out_path, root_label="Root")

    print(f"Reconstruction finished. Outputs in {cfg['output_dir']}/")
    return ancestors

def main():
    run(CONFIG)

if __name__ == "__main__":
    main()

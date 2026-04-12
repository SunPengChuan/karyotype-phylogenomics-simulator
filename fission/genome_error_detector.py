#!/usr/bin/env python3
"""
Genome Error Detector
=====================
Detects two classes of genomic errors by comparing extant species karyotypes
against reconstructed root ancestors:

  1. Misassembly  — a foreign gene block inserted into a target chromosome.
     Three detection methods:
       DUP      : duplicate genes forming a contiguous second-occurrence block
       SANDWICH : minority root-chromosome genes flanked by majority genes,
                  with both junction pairs absent from ancestor adjacencies
       SIB_DOM  : sibling-branch dominance — catches misassemblies that pass
                  the ancestor filter by comparing gene counts across tree branches

  2. Fission — a root chromosome split into two or more disjoint pure fragments
     in an extant species (WGD duplicates are excluded via overlap check).

Usage (run from the experiment directory):
  python genome_error_detector.py [options]
"""

import argparse
import os
import re
from collections import defaultdict, Counter


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def parse_ancestors(filepath):
    """Return {node: {chr_id: [genes]}} preserving gene order."""
    data = defaultdict(dict)
    with open(filepath, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('Node'):
                continue
            parts = line.split('\t')
            if len(parts) < 4:
                continue
            node, chr_id, genes_str = parts[0], parts[1], parts[3]
            data[node][chr_id] = genes_str.split()
    return dict(data)


def parse_species_karyotypes(filepath):
    """Return {species: {chr_name: [genes]}}."""
    data = defaultdict(dict)
    with open(filepath, encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3 and parts[0] not in ('Outgroup', 'Species', '---'):
                data[parts[0]][parts[1]] = parts[2:]
    return dict(data)


def parse_tree(newick_path):
    """
    Parse a Newick file with ete3 and return tree topology info.

    Returns a dict with:
      children_map  : {node_name: [child_names]}
      pre_wgd_nodes : set of nodes whose name ends with '_preWGD'
      post_wgd_nodes: set of all nodes that are descendants of a _preWGD node
                      (i.e. the WGD-derived subtrees to exclude from SIB_DOM)
    """
    from ete3 import Tree
    with open(newick_path, encoding='utf-8') as f:
        newick = f.read().strip()
    t = Tree(newick, format=1)

    children_map = {
        node.name: [c.name for c in node.children if c.name]
        for node in t.traverse() if node.name
    }

    pre_wgd_nodes = {n for n in children_map if n.endswith('_preWGD')}

    post_wgd_nodes = set()
    def _mark_subtree(name):
        post_wgd_nodes.add(name)
        for child in children_map.get(name, []):
            _mark_subtree(child)
    for pre in pre_wgd_nodes:
        for child in children_map.get(pre, []):
            if not child.endswith('_preWGD'):
                _mark_subtree(child)

    return {
        'children_map'  : children_map,
        'pre_wgd_nodes' : pre_wgd_nodes,
        'post_wgd_nodes': post_wgd_nodes,
    }


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def gene_family(g):
    """Return the alphabetic prefix of a gene name (e.g. 'E' for 'E292')."""
    m = re.match(r'[A-Za-z]+', g)
    return m.group() if m else ''


# ---------------------------------------------------------------------------
# Misassembly detection (methods A, B)
# ---------------------------------------------------------------------------

def detect_misassemblies(root_karyo, species_karyo, species_name,
                         ancestor_gsets, ancestor_adj,
                         min_dup_contiguity=0.7):
    """
    Detect misassemblies in one species using methods DUP and SANDWICH.

    DUP   : detects a contiguous block of duplicate genes.
    SANDWICH: detects a foreign gene block flanked by the same root chromosome
              on both sides, with both junction pairs absent from all ancestor
              chromosome adjacencies (novel junctions = assembly artifact).

    Returns a list of event dicts.
    """
    root_g2c = {}
    root_gene_sets = {}
    for cid, genes in root_karyo.items():
        root_gene_sets[cid] = set(genes)
        for g in genes:
            root_g2c[g] = cid

    results = []

    for sp_cid, genes in species_karyo.items():
        total = len(genes)
        if total < 10:
            continue

        # -- Method A: DUP (duplicate-gene block) ----------------------------
        gene_counts = Counter(genes)
        dups = {g for g, c in gene_counts.items() if c > 1}
        if dups:
            second_pos = {}
            for i, g in enumerate(genes):
                if g in dups and g not in second_pos:
                    if genes[:i].count(g) >= 1:
                        second_pos[g] = i

            if second_pos:
                positions = sorted(second_pos.values())
                span = positions[-1] - positions[0] + 1
                contiguity = len(dups) / span if span > 0 else 0

                if contiguity >= min_dup_contiguity:
                    s, e = positions[0], positions[-1]
                    insert_fams = Counter(gene_family(g) for g in genes[s:e + 1])
                    breakpoints = []
                    if s > 0:
                        breakpoints.append((genes[s - 1], genes[s]))
                    if e < total - 1:
                        breakpoints.append((genes[e], genes[e + 1]))
                    results.append({
                        'species'    : species_name,
                        'chr'        : sp_cid,
                        'method'     : 'DUP',
                        'insert_pos' : s,
                        'insert_size': len(dups),
                        'insert_fams': dict(insert_fams),
                        'breakpoints': breakpoints,
                        'contiguity' : round(contiguity, 2),
                    })
                    continue  # skip SANDWICH for this chromosome

        # -- Method B: SANDWICH (foreign block flanked by main-family genes) --
        assignments = [root_g2c.get(g) for g in genes]
        root_cnts = Counter(rc for rc in assignments if rc)
        if len(root_cnts) < 2:
            continue

        main_chr  = max(root_cnts, key=root_cnts.get)
        main_cnt  = root_cnts[main_chr]
        main_size = len(root_gene_sets.get(main_chr, []))

        # Filter 1: WGD — main count must not exceed root chromosome size
        if main_cnt > main_size:
            continue

        # Filter 2: ancestor — gene set already explained by reconstruction
        if frozenset(genes) in ancestor_gsets:
            continue

        # Find every maximal contiguous block that contains no main_chr genes.
        # A block is sandwiched when the same root chromosome flanks it on both
        # sides AND both junction pairs are absent from all ancestor adjacencies
        # (novel junctions indicate an assembly artifact, not a natural event).
        i = 0
        while i < total:
            if assignments[i] is not None and assignments[i] != main_chr:
                # Extend until we hit a main_chr gene or end of chromosome
                j = i
                while j < total and assignments[j] != main_chr:
                    j += 1
                bstart, bend = i, j - 1

                flanking_before = assignments[bstart - 1] if bstart > 0 else None
                flanking_after  = assignments[bend + 1]  if bend < total - 1 else None

                if (flanking_before and flanking_after
                        and flanking_before == flanking_after):
                    j_before = (genes[bstart - 1], genes[bstart])
                    j_after  = (genes[bend],        genes[bend + 1])

                    if (j_before not in ancestor_adj
                            and j_after not in ancestor_adj):
                        block_gset = frozenset(genes[bstart:bend + 1])
                        if block_gset not in ancestor_gsets:
                            insert_fams = Counter(
                                gene_family(genes[k])
                                for k in range(bstart, bend + 1)
                                if assignments[k]
                            )
                            minor_chrs = {
                                assignments[k]
                                for k in range(bstart, bend + 1)
                                if assignments[k]
                            }
                            results.append({
                                'species'    : species_name,
                                'chr'        : sp_cid,
                                'method'     : 'SANDWICH',
                                'insert_pos' : bstart,
                                'insert_size': bend - bstart + 1,
                                'insert_fams': dict(insert_fams),
                                'breakpoints': [j_before, j_after],
                                'main_chr'   : flanking_before,
                                'minor_chr'  : str(minor_chrs),
                                'ratio'      : round((bend - bstart + 1) / total * 100, 1),
                            })
                i = j if j > i else i + 1
            else:
                i += 1

    return results


# ---------------------------------------------------------------------------
# Sibling-branch dominance detection (method C)
# ---------------------------------------------------------------------------

def detect_sibling_dominance(root_karyo, species_data, ancestor_data,
                              ancestor_gsets, tree_info,
                              min_ratio=1.0, min_excess=5):
    """
    Detect misassemblies that passed the ancestor filter (Filter 2) by checking
    whether a root chromosome's gene count in the candidate chromosome exceeds
    the maximum found in any non-WGD sibling branch.

    A ratio > 1.0 with excess >= min_excess genes flags the chromosome as
    suspicious (the inflated count is unexplained by normal inheritance).

    Returns a list of event dicts.
    """
    children_map   = tree_info['children_map']
    pre_wgd_nodes  = tree_info['pre_wgd_nodes']
    post_wgd_nodes = tree_info['post_wgd_nodes']

    root_g2c = {g: cid for cid, genes in root_karyo.items() for g in genes}

    results = []

    for sp_name, sp_karyo in species_data.items():
        if sp_name == 'Outgroup':
            continue
        for sp_cid, genes in sp_karyo.items():
            if len(genes) < 10:
                continue

            root_cnts = Counter(rc for rc in (root_g2c.get(g) for g in genes) if rc)
            if len(root_cnts) < 2:
                continue

            # Only examine chromosomes that passed the ancestor filter
            gset = frozenset(genes)
            if gset not in ancestor_gsets:
                continue

            anc_node, anc_chr = ancestor_gsets[gset]
            if anc_node in post_wgd_nodes or anc_node in pre_wgd_nodes:
                continue

            # Compare against each non-WGD child branch of the ancestor node
            non_wgd_children = [
                c for c in children_map.get(anc_node, [])
                if c not in post_wgd_nodes and c not in pre_wgd_nodes
            ]
            if not non_wgd_children:
                continue

            flagged = []
            for child in non_wgd_children:
                child_data = ancestor_data.get(child) or species_data.get(child)
                if not child_data:
                    continue
                for rc, count_cand in root_cnts.items():
                    max_branch = max(
                        (sum(1 for g in cg if root_g2c.get(g) == rc)
                         for cg in child_data.values()),
                        default=0
                    )
                    if max_branch == 0:
                        continue
                    ratio  = count_cand / max_branch
                    excess = count_cand - max_branch
                    if ratio > min_ratio and excess >= min_excess:
                        flagged.append({
                            'root_chr'     : rc,
                            'count_in_cand': count_cand,
                            'max_in_branch': max_branch,
                            'branch'       : child,
                            'ratio'        : round(ratio, 2),
                            'excess'       : excess,
                        })

            if flagged:
                best = max(flagged, key=lambda x: x['ratio'])
                results.append({
                    'species'         : sp_name,
                    'chr'             : sp_cid,
                    'method'          : 'SIB_DOM',
                    'anc_node'        : anc_node,
                    'anc_chr'         : anc_chr,
                    'best_signal'     : best,
                    'root_composition': dict(root_cnts),
                })

    return results


# ---------------------------------------------------------------------------
# Fission detection
# ---------------------------------------------------------------------------

def detect_fissions(root_karyo, species_karyo, species_name):
    """
    Detect fission events for one species.

    A fission is confirmed when a root chromosome is covered by 2+ disjoint
    pure fragments (single-gene-family chromosomes) in the species.
    Overlapping fragments indicate WGD duplication and are excluded.

    Returns a list of event dicts.
    """
    root_gene_sets = {cid: set(genes) for cid, genes in root_karyo.items()}

    results = []

    for root_cid, root_genes in root_karyo.items():
        root_set   = root_gene_sets[root_cid]
        root_total = len(root_genes)

        pure_fragments = []
        complete_copies = []

        for sp_cid, sp_genes in species_karyo.items():
            overlap = set(sp_genes) & root_set
            if not overlap:
                continue
            if len(overlap) == root_total:
                complete_copies.append(sp_cid)
            elif len({gene_family(g) for g in sp_genes}) == 1:
                pure_fragments.append((sp_cid, overlap))

        if len(pure_fragments) < 2:
            continue

        # Fission fragments must be disjoint; overlap means WGD duplication
        frag_sets = [ov for _, ov in pure_fragments]
        if any(frag_sets[i] & frag_sets[j]
               for i in range(len(frag_sets))
               for j in range(i + 1, len(frag_sets))):
            continue

        covered = set().union(*frag_sets)
        coverage = len(covered) / root_total
        if coverage < 0.5:
            continue

        # Locate breakpoints in root gene order
        gene_to_frag = {g: idx for idx, (_, ov) in enumerate(pure_fragments)
                        for g in ov}
        breakpoints = []
        prev_g = prev_frag = None
        for g in root_genes:
            curr_frag = gene_to_frag.get(g)
            if curr_frag is None:
                prev_g, prev_frag = g, None
                continue
            if prev_frag is not None and curr_frag != prev_frag:
                breakpoints.append((prev_g, g))
            prev_g, prev_frag = g, curr_frag

        if not breakpoints:
            continue

        results.append({
            'species'       : species_name,
            'root_chr'      : root_cid,
            'root_size'     : root_total,
            'breakpoints'   : breakpoints,
            'pure_fragments': [(sp_cid, len(ov)) for sp_cid, ov in pure_fragments],
            'shared_copies' : len(complete_copies),
            'coverage_pct'  : round(coverage * 100, 1),
        })

    return results


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def generate_report(misassemblies, fissions, sib_dom, output_path):
    """Write a human-readable detection report to output_path."""
    lines = []
    sep = '=' * 70

    lines += [sep, 'GENOME ERROR DETECTOR', sep,
              f'Misassemblies detected : {len(misassemblies) + len(sib_dom)}',
              f'Fission events detected: {len(fissions)}', '']

    if misassemblies:
        lines += [sep, 'MISASSEMBLIES (DUP / SANDWICH)', sep]
        for ev in misassemblies:
            lines.append(f"Species: {ev['species']}  Chr: {ev['chr']}  Method: {ev['method']}")
            lines.append(f"  Insert position: {ev['insert_pos']}  Insert size: {ev['insert_size']} genes")
            lines.append(f"  Insert families: {ev['insert_fams']}")
            if ev['method'] == 'SANDWICH':
                lines.append(f"  Main root chr: {ev['main_chr']}  "
                             f"Insert root chr: {ev['minor_chr']} ({ev['ratio']}%)")
            else:
                lines.append(f"  Duplicate contiguity: {ev['contiguity']}")
            lines.append('  Breakpoints:')
            for gl, gr in ev['breakpoints']:
                lines.append(f"    {gl} | {gr}")
            lines.append('')

    if sib_dom:
        lines += [sep, 'MISASSEMBLIES (SIB_DOM)', sep]
        for ev in sib_dom:
            sig = ev['best_signal']
            lines.append(f"Species: {ev['species']}  Chr: {ev['chr']}  Method: SIB_DOM")
            lines.append(f"  Ancestor node: {ev['anc_node']}  Ancestor chr: {ev['anc_chr']}")
            lines.append(f"  Anomalous root chr: {sig['root_chr']}  "
                         f"Count in candidate: {sig['count_in_cand']}  "
                         f"Max in sibling branch: {sig['max_in_branch']}  "
                         f"Ratio: {sig['ratio']}x  Excess: {sig['excess']}")
            lines.append(f"  Root composition: {ev['root_composition']}")
            lines.append('')

    if fissions:
        lines += [sep, 'FISSION EVENTS', sep]
        for ev in fissions:
            lines.append(f"Species: {ev['species']}")
            lines.append(f"  Root chr   : {ev['root_chr']}  ({ev['root_size']} genes)")
            lines.append(f"  Full copies: {ev['shared_copies']}")
            lines.append('  Pure fragments:')
            for sp_cid, cnt in sorted(ev['pure_fragments']):
                lines.append(f"    {sp_cid}  ({cnt} genes)")
            lines.append('  Breakpoints (root gene order):')
            for gl, gr in ev['breakpoints']:
                lines.append(f"    {gl} | {gr}")
            lines.append(f"  Fragment coverage: {ev['coverage_pct']}%")
            lines.append('')

    # Summary table
    lines += [sep, 'SUMMARY TABLE', sep]
    if misassemblies:
        lines += ['--- Misassemblies ---',
                  f"{'Species':<8} {'Chr':<8} {'Method':<10} {'Pos':<8} {'Size':<8} Breakpoints",
                  '-' * 80]
        for ev in misassemblies:
            bp = '; '.join(f"{gl}|{gr}" for gl, gr in ev['breakpoints'])
            lines.append(f"{ev['species']:<8} {ev['chr']:<8} {ev['method']:<10} "
                         f"{ev['insert_pos']:<8} {ev['insert_size']:<8} {bp}")
    if sib_dom:
        lines += ['--- Misassemblies (SIB_DOM) ---',
                  f"{'Species':<8} {'Chr':<8} {'AncNode':<10} {'RootChr':<16} {'Ratio':<8} Excess",
                  '-' * 80]
        for ev in sib_dom:
            sig = ev['best_signal']
            lines.append(f"{ev['species']:<8} {ev['chr']:<8} {ev['anc_node']:<10} "
                         f"{sig['root_chr']:<16} {sig['ratio']:<8} {sig['excess']}")
    if fissions:
        lines += ['--- Fission events ---',
                  f"{'Species':<8} {'RootChr':<14} {'Breakpoints':<30} Fragments",
                  '-' * 80]
        for ev in fissions:
            bp   = '; '.join(f"{gl}|{gr}" for gl, gr in ev['breakpoints'])
            frgs = ', '.join(f"{c}({n})" for c, n in sorted(ev['pure_fragments']))
            lines.append(f"{ev['species']:<8} {ev['root_chr']:<14} {bp:<30} {frgs}")

    report = '\n'.join(lines)
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(report)
    return report


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Detect misassembly errors and fission events in karyotype data.')
    parser.add_argument('--ancestors', default='output_reconstruction_fission/ancestor_gene_sets_by_node.tsv',
        help='Reconstructed ancestor gene sets TSV (default: %(default)s)')
    parser.add_argument('--species', default='output_simulator/karyotypes_fission.txt',
        help='Extant species karyotype file (default: %(default)s)')
    parser.add_argument('--tree', default='output_simulator/tree.nwk',
        help='Newick tree file for SIB_DOM detection (default: %(default)s)')
    parser.add_argument('--output', default='error_detection_report.txt',
        help='Output report file (default: %(default)s)')
    parser.add_argument('--min-dup-contiguity', type=float, default=0.7,
        help='Min contiguity for DUP detection (default: %(default)s)')
    args = parser.parse_args()

    print('Loading data...')
    ancestor_data = parse_ancestors(args.ancestors)
    species_data  = parse_species_karyotypes(args.species)

    root_karyo = ancestor_data.get('Root')
    if not root_karyo:
        print('Error: Root node not found in ancestor TSV')
        return

    # Build ancestor gene-set index for filtering natural rearrangements
    ancestor_gsets = {}
    for node, node_data in ancestor_data.items():
        if node == 'Root':
            continue
        for cid, genes in node_data.items():
            ancestor_gsets[frozenset(genes)] = (node, cid)

    # Build ancestor adjacency set: all adjacent gene pairs across all ancestor chromosomes.
    # Novel junctions (absent from this set) are a strong signal for assembly artifacts.
    ancestor_adj = set()
    for node_data in ancestor_data.values():
        for genes in node_data.values():
            for i in range(len(genes) - 1):
                ancestor_adj.add((genes[i], genes[i + 1]))

    print(f'Root: {len(root_karyo)} chromosomes')
    print(f'Ancestor nodes: {len(ancestor_data) - 1} (excluding Root)')
    print()

    all_misassemblies = []
    all_fissions      = []

    for sp_name in sorted(species_data):
        if sp_name == 'Outgroup':
            continue
        sp_karyo = species_data[sp_name]

        for ev in detect_misassemblies(root_karyo, sp_karyo, sp_name, ancestor_gsets,
                                       ancestor_adj, args.min_dup_contiguity):
            all_misassemblies.append(ev)
            bp = '; '.join(f"{gl}|{gr}" for gl, gr in ev['breakpoints'])
            print(f"[MISASSEMBLY/{ev['method']}] {sp_name}:{ev['chr']}  "
                  f"pos={ev['insert_pos']}  size={ev['insert_size']}  "
                  f"fams={ev['insert_fams']}")
            print(f"  breakpoints: {bp}")

        for ev in detect_fissions(root_karyo, sp_karyo, sp_name):
            all_fissions.append(ev)
            bp   = '; '.join(f"{gl}|{gr}" for gl, gr in ev['breakpoints'])
            frgs = ', '.join(f"{c}({n})" for c, n in sorted(ev['pure_fragments']))
            print(f"[FISSION]      {sp_name}  root={ev['root_chr']}  breakpoints: {bp}")
            print(f"  fragments: {frgs}  coverage: {ev['coverage_pct']}%")

    # Method C: sibling-branch dominance (requires tree file)
    sib_dom = []
    if os.path.exists(args.tree):
        tree_info = parse_tree(args.tree)
        sib_dom = detect_sibling_dominance(
            root_karyo, species_data, ancestor_data, ancestor_gsets, tree_info)
        for ev in sib_dom:
            sig = ev['best_signal']
            print(f"[MISASSEMBLY/SIB_DOM] {ev['species']}:{ev['chr']}  "
                  f"anc={ev['anc_node']}/{ev['anc_chr']}  "
                  f"root_chr={sig['root_chr']}  "
                  f"count={sig['count_in_cand']}  max_branch={sig['max_in_branch']}  "
                  f"ratio={sig['ratio']}x  excess={sig['excess']}")
    else:
        print(f'Warning: tree file not found ({args.tree}), skipping SIB_DOM detection')

    total_ma = len(all_misassemblies) + len(sib_dom)
    print(f'\nDetected {total_ma} misassembl{"y" if total_ma == 1 else "ies"} '
          f'({len(sib_dom)} via SIB_DOM), '
          f'{len(all_fissions)} fission event{"" if len(all_fissions) == 1 else "s"}')
    generate_report(all_misassemblies, all_fissions, sib_dom, args.output)
    print(f'Report saved: {args.output}')


if __name__ == '__main__':
    main()

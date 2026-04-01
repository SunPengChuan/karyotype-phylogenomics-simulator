#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import random
import subprocess
import sys
from collections import Counter
import shutil

# Source file paths (relative to this script)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SIMULATOR_SRC = os.path.join(SCRIPT_DIR, "evolution_simulator.py")
RECONSTRUCTION_SRC = os.path.join(SCRIPT_DIR, "ancestor_reconstruction.py")

PROB_RANGES = {
    "inversion_prob": (0.5, 0.95),
    "translocation_rct": (0.05, 0.8),
    "fusion_ncf": (0.0, 0.6),
    "fusion_eej": (0.0, 0.8),
    "wgd_probability": (0.0, 0.2),
}


def _parse_int_list(value: str):
    value = value.strip()
    if not value:
        return []
    parts = [p.strip() for p in value.split(",") if p.strip()]
    out = []
    for p in parts:
        if ".." in p:
            a, b = p.split("..", 1)
            out.extend(list(range(int(a), int(b) + 1)))
        else:
            out.append(int(p))
    return out


def _ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)
    return path


def _remove_path(path: str):
    if not os.path.exists(path):
        return
    if os.path.isdir(path):
        shutil.rmtree(path)
    else:
        os.remove(path)


import re


def _modify_simulator_config(src_path: str, dst_path: str, sim_params: dict):
    """Modify CONFIG parameters in evolution_simulator.py"""
    with open(src_path, "r", encoding="utf-8") as f:
        content = f.read()
    
    # Use regex to replace parameters
    def replace_param(content: str, param_name: str, new_value) -> str:
        """Replace parameter value using regex"""
        if isinstance(new_value, str):
            pattern = rf'({param_name}\s*=\s*")[^"]*(")'
            # Use lambda to prevent backslashes in paths from being treated as escape characters
            replacement = lambda m: m.group(1) + new_value.replace('\\', '\\\\') + m.group(2)
            return re.sub(pattern, replacement, content)
        elif isinstance(new_value, bool):
            pattern = rf'({param_name}\s*=\s*)(True|False)'
            replacement = rf'\g<1>{new_value}'
        elif isinstance(new_value, (int, float)):
            pattern = rf'({param_name}\s*=\s*)[\d.]+'
            replacement = rf'\g<1>{new_value}'
        else:
            return content
        return re.sub(pattern, replacement, content)
    
    # Replace simulator parameters
    content = replace_param(content, "seed", sim_params["seed"])
    content = replace_param(content, "num_modern_species", sim_params["num_modern_species"])
    content = replace_param(content, "num_ancestor_chromosomes", sim_params["num_ancestor_chromosomes"])
    content = replace_param(content, "save_dir", sim_params["save_dir"])
    content = replace_param(content, "visualize_collinearity", sim_params.get("visualize_collinearity", False))
    if sim_params.get("inversion_prob") is not None:
        content = replace_param(content, "inversion_prob", float(sim_params["inversion_prob"]))
    if sim_params.get("translocation_rct") is not None:
        content = replace_param(content, "translocation_rct", float(sim_params["translocation_rct"]))
    if sim_params.get("fusion_ncf") is not None:
        content = replace_param(content, "fusion_ncf", float(sim_params["fusion_ncf"]))
    if sim_params.get("fusion_eej") is not None:
        content = replace_param(content, "fusion_eej", float(sim_params["fusion_eej"]))
    if sim_params.get("wgd_probability") is not None:
        content = replace_param(content, "wgd_probability", float(sim_params["wgd_probability"]))
    
    with open(dst_path, "w", encoding="utf-8") as f:
        f.write(content)


def _parse_events(events_path: str):
    counts = Counter()
    if not os.path.exists(events_path):
        return counts
    with open(events_path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    for line in lines[2:]:
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        etype = parts[1].strip()
        if etype == "PreWGD":
            continue
        counts[etype] += 1
    return counts


def _format_count_and_rate(count: int, branch_count: int):
    rate = 0.0 if branch_count <= 0 else (count / branch_count)
    return f"{count} ({rate:.2f})"


def _count_branches_from_nwk(nwk_path: str) -> int:
    """Count branches from Newick file"""
    if not os.path.exists(nwk_path):
        return 0
    try:
        from ete3 import Tree
        t = Tree(nwk_path, format=1)
        return sum(1 for n in t.traverse("preorder") if not n.is_root())
    except Exception:
        return 0


def _count_species_from_nwk(nwk_path: str) -> int:
    """Count species from Newick file"""
    if not os.path.exists(nwk_path):
        return 0
    try:
        from ete3 import Tree
        t = Tree(nwk_path, format=1)
        return len(t.get_leaf_names())
    except Exception:
        return 0


def _load_true_root_karyotype(true_root_path: str) -> dict:
    """Load true root karyotype from file"""
    if not true_root_path or not os.path.exists(true_root_path):
        return {}
    
    karyo = {}
    with open(true_root_path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    
    for line in lines:
        parts = line.strip().split()
        if len(parts) < 3:
            continue
        sp_name = parts[0]
        if sp_name != "TrueRoot":
            continue
        cid = parts[1]
        genes = parts[2:]
        karyo[cid] = set(genes)
    
    return karyo


def _load_reconstructed_root(reconstructed_tsv_path: str) -> dict:
    """Load reconstructed root karyotype from TSV file"""
    if not reconstructed_tsv_path or not os.path.exists(reconstructed_tsv_path):
        return {}
    
    karyo = {}
    with open(reconstructed_tsv_path, "r", encoding="utf-8") as f:
        header_skipped = False
        for line in f:
            line = line.strip()
            if not line:
                continue
            if not header_skipped:
                header_skipped = True
                continue
            parts = line.split("\t", 3)
            if len(parts) < 4:
                continue
            node, cid, _, genes_str = parts
            if node != "Root":
                continue
            genes = genes_str.split()
            karyo[cid] = set(genes)
    
    return karyo


def _calculate_reconstruction_accuracy(true_root_path: str, reconstructed_tsv_path: str) -> dict:
    """Calculate reconstruction accuracy by comparing true root and reconstructed root.
    
    Returns a dict with:
        - true_chr_count: number of true chromosomes
        - reconstructed_chr_count: number of reconstructed chromosomes
        - exact_match_count: number of true chromosomes exactly matched
        - exact_match_rate: exact_match_count / true_chr_count
        - gene_coverage: fraction of true root genes covered by reconstructed chromosomes
    """
    true_root = _load_true_root_karyotype(true_root_path)
    rec_root = _load_reconstructed_root(reconstructed_tsv_path)
    
    if not true_root:
        return {"true_chr_count": 0, "reconstructed_chr_count": 0, 
                "exact_match_count": 0, "exact_match_rate": 0.0, "gene_coverage": 0.0}
    
    if not rec_root:
        return {"true_chr_count": len(true_root), "reconstructed_chr_count": 0,
                "exact_match_count": 0, "exact_match_rate": 0.0, "gene_coverage": 0.0}
    
    # Calculate exact matches
    exact_matches = 0
    for cid, true_genes in true_root.items():
        for rec_cid, rec_genes in rec_root.items():
            if true_genes == rec_genes:
                exact_matches += 1
                break
    
    # Calculate gene coverage
    all_true_genes = set()
    for genes in true_root.values():
        all_true_genes.update(genes)
    
    all_rec_genes = set()
    for genes in rec_root.values():
        all_rec_genes.update(genes)
    
    gene_coverage = len(all_true_genes & all_rec_genes) / len(all_true_genes) if all_true_genes else 0.0
    
    return {
        "true_chr_count": len(true_root),
        "reconstructed_chr_count": len(rec_root),
        "exact_match_count": exact_matches,
        "exact_match_rate": exact_matches / len(true_root) if true_root else 0.0,
        "gene_coverage": gene_coverage
    }


def _write_tsv(path: str, header, rows):
    _ensure_dir(os.path.dirname(path))
    with open(path, "w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for row in rows:
            f.write("\t".join(str(row.get(col, "")) for col in header) + "\n")


def _extract_simulator_probabilities(sim_path: str) -> dict:
    out = {
        "inversion_prob": "",
        "translocation_rct": "",
        "fusion_ncf": "",
        "fusion_eej": "",
        "wgd_probability": "",
    }
    if not sim_path or not os.path.exists(sim_path):
        return out
    try:
        with open(sim_path, "r", encoding="utf-8") as f:
            content = f.read()
    except OSError:
        return out

    def _grab(name: str):
        m = re.search(rf"\b{name}\s*=\s*([0-9]*\.?[0-9]+)", content)
        return "" if not m else m.group(1)

    for k in list(out.keys()):
        out[k] = _grab(k)
    return out


def _scenario_is_done(scenario_dir: str) -> bool:
    sim_dir = os.path.join(scenario_dir, "output_simulator")
    recon_dir = os.path.join(scenario_dir, "output_reconstruction")
    events_path = os.path.join(sim_dir, "events.txt")
    tree_nwk_path = os.path.join(sim_dir, "tree.nwk")
    true_root_path = os.path.join(sim_dir, "karyotypes_true_root.txt")
    reconstructed_tsv_path = os.path.join(recon_dir, "ancestor_gene_sets_by_node.tsv")
    return all(
        os.path.exists(p)
        for p in (events_path, tree_nwk_path, true_root_path, reconstructed_tsv_path)
    )


def run_scenarios(
    num_scenarios,
    species_values,
    anc_chr_values,
    out_dir,
    param_mode="cycle",
    base_seed=100,
    skip_run=False,
    skip_existing=False,
    visualize_collinearity=True,
):
    """Run multiple test scenarios, each with independent program copies
    
    Args:
        skip_run: If True, only collect results from existing directories without running simulations
        skip_existing: If True, skip running scenarios whose outputs already exist
    """
    results = []

    species_seq = list(species_values)
    anc_chr_seq = list(anc_chr_values)
    if param_mode == "shuffle":
        rng = random.Random(base_seed)
        rng.shuffle(species_seq)
        rng.shuffle(anc_chr_seq)

    def _rand_prob(rng: random.Random, key: str) -> float:
        lo, hi = PROB_RANGES[key]
        v = rng.uniform(float(lo), float(hi))
        if v < 0.0:
            v = 0.0
        if v > 1.0:
            v = 1.0
        return float(f"{v:.3f}")

    for idx in range(num_scenarios):
        seed = base_seed + idx
        if param_mode == "random":
            rng = random.Random(seed)
            requested_species = rng.choice(species_seq)
            num_anc_chr = rng.choice(anc_chr_seq)
        else:
            requested_species = species_seq[idx % len(species_seq)]
            num_anc_chr = anc_chr_seq[(idx // len(species_seq)) % len(anc_chr_seq)]

        scenario_dir = os.path.join(out_dir, f"scenario_{idx+1:03d}")
        sim_dir = os.path.join(scenario_dir, "output_simulator")
        recon_dir = os.path.join(scenario_dir, "output_reconstruction")
        reconstructed_tsv_path = os.path.join(recon_dir, "ancestor_gene_sets_by_node.tsv")
        true_root_path = os.path.join(sim_dir, "karyotypes_true_root.txt")
        tree_nwk_path = os.path.join(sim_dir, "tree.nwk")
        events_path = os.path.join(sim_dir, "events.txt")
        sim_path = os.path.join(scenario_dir, "evolution_simulator.py")

        done = _scenario_is_done(scenario_dir)
        should_run = (not skip_run) and (not (skip_existing and done))

        if should_run:
            _ensure_dir(scenario_dir)
            _ensure_dir(sim_dir)

            # 1. Copy source files to scenario directory
            sim_dst = os.path.join(scenario_dir, "evolution_simulator.py")
            recon_dst = os.path.join(scenario_dir, "ancestor_reconstruction.py")
            shutil.copy(SIMULATOR_SRC, sim_dst)
            shutil.copy(RECONSTRUCTION_SRC, recon_dst)

            rng_p = random.Random(seed + 100003)
            p_inv = _rand_prob(rng_p, "inversion_prob")
            p_rct = _rand_prob(rng_p, "translocation_rct")
            p_ncf = _rand_prob(rng_p, "fusion_ncf")
            p_eej = _rand_prob(rng_p, "fusion_eej")
            p_wgd = _rand_prob(rng_p, "wgd_probability")

            # 2. Modify simulator config
            sim_params = {
                "seed": seed,
                "num_modern_species": requested_species,
                "num_ancestor_chromosomes": num_anc_chr,
                "save_dir": "output_simulator",
                "visualize_collinearity": visualize_collinearity,
                "inversion_prob": p_inv,
                "translocation_rct": p_rct,
                "fusion_ncf": p_ncf,
                "fusion_eej": p_eej,
                "wgd_probability": p_wgd,
            }
            _modify_simulator_config(sim_dst, sim_dst, sim_params)

        # 3. Run tests (skip if skip_run is True)
        print(f"\n{'='*60}")
        print(f"Scenario {idx+1}/{num_scenarios}: species={requested_species}, anc_chr={num_anc_chr}, seed={seed}")
        print(f"Directory: {scenario_dir}")
        print(f"{'='*60}")
        
        if skip_run:
            print("SKIP: Skipping simulation and reconstruction (skip_run=True)")
        elif skip_existing and done:
            print("SKIP: Outputs already exist for this scenario (skip_existing=True)")
        else:
            try:
                result = subprocess.run(
                    [sys.executable, "evolution_simulator.py"],
                    cwd=scenario_dir,
                    timeout=300,
                )
                if result.returncode != 0:
                    print(f"ERROR: Simulator returned code {result.returncode}")
                    continue
            except subprocess.TimeoutExpired:
                print(f"TIMEOUT: Scenario {idx+1} simulator exceeded 5 minutes")
                continue
            except Exception as e:
                print(f"ERROR running scenario {idx+1} simulator: {e}")
                continue

            # 4. Run ancestor_reconstruction.py
            try:
                result = subprocess.run(
                    [sys.executable, "ancestor_reconstruction.py"],
                    cwd=scenario_dir,
                    timeout=300,
                )
                if result.returncode != 0:
                    print(f"ERROR: Reconstruction returned code {result.returncode}")
            except subprocess.TimeoutExpired:
                print(f"TIMEOUT: Scenario {idx+1} reconstruction exceeded 5 minutes")
            except Exception as e:
                print(f"ERROR running scenario {idx+1} reconstruction: {e}")

        # 5. Collect results
        event_counts = _parse_events(events_path)
        
        branch_count = _count_branches_from_nwk(tree_nwk_path)
        num_species = _count_species_from_nwk(tree_nwk_path)
        
        # Calculate reconstruction accuracy
        accuracy = _calculate_reconstruction_accuracy(true_root_path, reconstructed_tsv_path)

        sim_probs = _extract_simulator_probabilities(sim_path)
        
        # Print accuracy info
        if accuracy["true_chr_count"] > 0:
            print(f"Reconstruction Accuracy: {accuracy['exact_match_count']}/{accuracy['true_chr_count']} "
                  f"({accuracy['exact_match_rate']*100:.1f}%) chromosomes matched, "
                  f"gene coverage: {accuracy['gene_coverage']*100:.1f}%")

        results.append(
            dict(
                Scenario=f"scenario_{idx+1:03d}",
                SpeciesCount=num_species,
                AncestorChromosomeCount=num_anc_chr,
                BranchCount=branch_count,
                TrueChrCount=accuracy["true_chr_count"],
                ReconstructedChrCount=accuracy["reconstructed_chr_count"],
                ExactMatchCount=accuracy["exact_match_count"],
                ExactMatchRate=f"{accuracy['exact_match_rate']*100:.1f}%",
                GeneCoverage=f"{accuracy['gene_coverage']*100:.1f}%",
                P_Inversion=sim_probs.get("inversion_prob", ""),
                P_RCT=sim_probs.get("translocation_rct", ""),
                P_NCF=sim_probs.get("fusion_ncf", ""),
                P_EEJ=sim_probs.get("fusion_eej", ""),
                P_WGD=sim_probs.get("wgd_probability", ""),
                Inversion=_format_count_and_rate(event_counts.get("Inversion", 0), branch_count),
                RCT=_format_count_and_rate(event_counts.get("RCT", 0), branch_count),
                NCF=_format_count_and_rate(event_counts.get("NCF", 0), branch_count),
                EEJ=_format_count_and_rate(event_counts.get("EEJ", 0), branch_count),
                WGD=_format_count_and_rate(event_counts.get("WGD", 0), branch_count),
            )
        )

    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--num-scenarios", type=int, default=15)
    parser.add_argument("--species", type=str, default="8..20")
    parser.add_argument("--anc-chr", type=str, default="7..30")
    parser.add_argument("--out", type=str, default="output_experiments")
    parser.add_argument("--param-mode", type=str, default="random", choices=["cycle", "shuffle", "random"])
    parser.add_argument("--base-seed", type=int, default=42, help="Base random seed")
    parser.add_argument("--skip-run", action="store_true", help="Only collect results, skip running simulations")
    parser.add_argument("--skip-existing", action="store_true", help="Skip scenarios whose outputs already exist")
    parser.add_argument("--no-visualize", action="store_true", help="Disable collinearity visualization")
    args = parser.parse_args()

    species_values = _parse_int_list(args.species)
    anc_chr_values = _parse_int_list(args.anc_chr)
    if not species_values or not anc_chr_values:
        raise SystemExit("species and anc-chr must be non-empty")

    results = run_scenarios(
        num_scenarios=args.num_scenarios,
        species_values=species_values,
        anc_chr_values=anc_chr_values,
        out_dir=args.out,
        param_mode=args.param_mode,
        base_seed=args.base_seed,
        skip_run=args.skip_run,
        skip_existing=args.skip_existing,
        visualize_collinearity=not args.no_visualize,
    )

    header = [
        "Scenario",
        "SpeciesCount",
        "AncestorChromosomeCount",
        "BranchCount",
        "TrueChrCount",
        "ReconstructedChrCount",
        "ExactMatchCount",
        "ExactMatchRate",
        "GeneCoverage",
        "P_Inversion",
        "P_RCT",
        "P_NCF",
        "P_EEJ",
        "P_WGD",
        "Inversion",
        "RCT",
        "NCF",
        "EEJ",
        "WGD",
    ]

    _write_tsv(os.path.join(args.out, "results.tsv"), header, results)
    print(f"\nResults saved to: {os.path.join(args.out, 'results.tsv')}")


if __name__ == "__main__":
    main()

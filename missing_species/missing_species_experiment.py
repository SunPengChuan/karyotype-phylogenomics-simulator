#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch robustness test: simulate datasets with random species counts (13-23),
drop a random number of species (1 or 2), reconstruct, and report accuracy.

Usage:
    python missing_species_experiment.py              # 10 scenarios
    python missing_species_experiment.py --n 5
    python missing_species_experiment.py --skip-existing
"""

import argparse
import os
import random
import re
import shutil
import subprocess
import sys
from pathlib import Path

from ete3 import Tree

SCRIPT_DIR = Path(__file__).parent
SIMULATOR_SRC = SCRIPT_DIR / "evolution_simulator.py"
RECON_SRC     = SCRIPT_DIR / "ancestor_reconstruction.py"
DROP_SRC      = SCRIPT_DIR / "drop_species.py"

SPECIES_RANGE = (13, 23)   # inclusive


# ── helpers ───────────────────────────────────────────────────────────────────

def _ensure(p):
    Path(p).mkdir(parents=True, exist_ok=True)


def _count_leaves(nwk):
    if not os.path.exists(nwk):
        return 0
    try:
        return len(Tree(nwk, format=1).get_leaf_names())
    except Exception:
        return 0


def _load_true_root(path):
    karyo = {}
    if not os.path.exists(path):
        return karyo
    with open(path, encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3 and parts[0] == "TrueRoot":
                karyo[parts[1]] = set(parts[2:])
    return karyo


def _load_recon_root(path):
    karyo = {}
    if not os.path.exists(path):
        return karyo
    with open(path, encoding="utf-8") as f:
        next(f, None)
        for line in f:
            parts = line.strip().split("\t", 3)
            if len(parts) == 4 and parts[0] == "Root":
                karyo[parts[1]] = set(parts[3].split())
    return karyo


def _accuracy(true_root_path, recon_tsv_path):
    true  = _load_true_root(true_root_path)
    recon = _load_recon_root(recon_tsv_path)
    if not true:
        return dict(true_n=0, recon_n=0, exact=0, exact_rate=0.0, coverage=0.0)
    exact = sum(1 for genes in true.values() if genes in recon.values())
    all_true  = set().union(*true.values())
    all_recon = set().union(*recon.values()) if recon else set()
    return dict(
        true_n=len(true),
        recon_n=len(recon),
        exact=exact,
        exact_rate=exact / len(true),
        coverage=len(all_true & all_recon) / len(all_true) if all_true else 0.0,
    )


def _modify_sim_config(src, dst, params):
    content = Path(src).read_text(encoding="utf-8")
    for k, v in params.items():
        if isinstance(v, bool):
            content = re.sub(rf"({k}\s*=\s*)(True|False)", rf"\g<1>{v}", content)
        elif isinstance(v, str):
            content = re.sub(rf'({k}\s*=\s*")[^"]*(")', rf'\g<1>{v}\g<2>', content)
        else:
            content = re.sub(rf"({k}\s*=\s*)[\d.]+", rf"\g<1>{v}", content)
    Path(dst).write_text(content, encoding="utf-8")


def _modify_recon_config(src, dst, sim_dir, recon_dir):
    content = Path(src).read_text(encoding="utf-8")
    content = re.sub(r'input_tree\s*=\s*"[^"]*"',
                     f'input_tree="{sim_dir}/tree.nwk"', content)
    content = re.sub(r'input_karyotypes\s*=\s*"[^"]*"',
                     f'input_karyotypes="{sim_dir}/karyotypes_species_with_outgroup.txt"', content)
    content = re.sub(r'input_true_root_karyotype\s*=\s*"[^"]*"',
                     f'input_true_root_karyotype="{sim_dir}/karyotypes_true_root.txt"', content)
    content = re.sub(r'output_dir\s*=\s*"[^"]*"',
                     f'output_dir="{recon_dir}"', content)
    content = re.sub(r'tree_viz_output_dir\s*=\s*"[^"]*"',
                     f'tree_viz_output_dir="{sim_dir}"', content)
    Path(dst).write_text(content, encoding="utf-8")


def _get_species_list(karyotype_path):
    """Return base species names (no _preWGD, no Outgroup) from tab/space table."""
    seen, species = set(), []
    with open(karyotype_path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("Species") or line.startswith("-"):
                continue
            sp = line.split()[0]
            if sp not in ("Outgroup",) and "_preWGD" not in sp and sp not in seen:
                seen.add(sp)
                species.append(sp)
    return species


# ── main loop ─────────────────────────────────────────────────────────────────

def run_experiments(n, out_dir, base_seed, skip_existing):
    out_dir = Path(out_dir).resolve()
    rows = []

    for idx in range(n):
        seed = base_seed + idx
        rng  = random.Random(seed)

        # randomise species count and drop count per scenario
        num_species = rng.randint(*SPECIES_RANGE)
        n_drop      = rng.randint(1, 2)

        scenario_dir      = out_dir / f"scenario_{idx+1:03d}"
        sim_dir           = scenario_dir / "output_simulator"
        sim_dir_reduced   = scenario_dir / "output_simulator_reduced"
        recon_dir_reduced = scenario_dir / "output_reconstruction_reduced"
        recon_tsv         = recon_dir_reduced / "ancestor_gene_sets_by_node.tsv"
        true_root_path    = sim_dir / "karyotypes_true_root.txt"

        done = recon_tsv.exists()
        if skip_existing and done:
            print(f"\n[{idx+1}/{n}] Scenario {idx+1}: skipping (already done)")
        else:
            print(f"\n{'='*60}")
            print(f"[{idx+1}/{n}] seed={seed}  species={num_species}  drop={n_drop}")
            print(f"{'='*60}")

            _ensure(scenario_dir)

            # 1. simulate
            sim_script = scenario_dir / "evolution_simulator.py"
            shutil.copy(SIMULATOR_SRC, sim_script)
            _modify_sim_config(sim_script, sim_script, {
                "seed": seed,
                "num_modern_species": num_species,
                "save_dir": "output_simulator",
                "visualize_collinearity": False,
            })
            subprocess.run([sys.executable, "evolution_simulator.py"],
                           cwd=scenario_dir, check=True, timeout=300)

            # 2. drop species
            karyo_path   = sim_dir / "karyotypes_species_with_outgroup.txt"
            species_list = _get_species_list(str(karyo_path))
            dropped      = rng.sample(species_list, min(n_drop, len(species_list)))
            print(f"  Dropping: {dropped}")

            subprocess.run(
                [sys.executable, str(DROP_SRC),
                 "--drop", *dropped,
                 "--input-dir",  str(sim_dir),
                 "--output-dir", str(sim_dir_reduced)],
                cwd=scenario_dir, check=True, timeout=60,
            )

            # 3. reconstruct reduced only
            recon_script = scenario_dir / "recon_reduced.py"
            _modify_recon_config(RECON_SRC, recon_script,
                                 "output_simulator_reduced",
                                 "output_reconstruction_reduced")
            subprocess.run([sys.executable, "recon_reduced.py"],
                           cwd=scenario_dir, check=True, timeout=300)

            (scenario_dir / "dropped_species.txt").write_text(
                "\n".join(dropped) + "\n", encoding="utf-8")

        # collect results
        dropped_saved = []
        dp = scenario_dir / "dropped_species.txt"
        if dp.exists():
            dropped_saved = dp.read_text(encoding="utf-8").strip().splitlines()

        acc       = _accuracy(str(true_root_path), str(recon_tsv))
        n_full    = _count_leaves(str(sim_dir         / "tree.nwk"))
        n_reduced = _count_leaves(str(sim_dir_reduced / "tree.nwk"))

        print(f"  Species full/reduced: {n_full}/{n_reduced}  "
              f"Exact: {acc['exact']}/{acc['true_n']}  "
              f"Coverage: {acc['coverage']*100:.1f}%")

        rows.append(dict(
            Scenario=f"scenario_{idx+1:03d}",
            Seed=seed,
            SpeciesFull=n_full,
            SpeciesReduced=n_reduced,
            DroppedCount=len(dropped_saved),
            DroppedSpecies=",".join(dropped_saved),
            TrueChrCount=acc["true_n"],
            ReconChrCount=acc["recon_n"],
            ExactMatch=acc["exact"],
            ExactMatchRate=f"{acc['exact_rate']*100:.1f}%",
            GeneCoverage=f"{acc['coverage']*100:.1f}%",
        ))

    # write summary
    header = ["Scenario","Seed","SpeciesFull","SpeciesReduced","DroppedCount",
              "DroppedSpecies","TrueChrCount","ReconChrCount",
              "ExactMatch","ExactMatchRate","GeneCoverage"]
    tsv_path = out_dir / "summary.tsv"
    with open(tsv_path, "w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in header) + "\n")

    if rows:
        rates = [float(r["ExactMatchRate"].rstrip("%")) for r in rows]
        covs  = [float(r["GeneCoverage"].rstrip("%"))   for r in rows]
        print(f"\nSummary ({len(rows)} scenarios):")
        print(f"  Avg exact match rate : {sum(rates)/len(rates):.1f}%")
        print(f"  Avg gene coverage    : {sum(covs)/len(covs):.1f}%")
    print(f"Results saved to {tsv_path}")
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n",    type=int, default=10,
                        help="Number of scenarios (default 10)")
    parser.add_argument("--out",  default="output_experiments_missing",
                        help="Output directory")
    parser.add_argument("--seed", type=int, default=42,
                        help="Base random seed")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip scenarios whose outputs already exist")
    args = parser.parse_args()

    run_experiments(
        n=args.n,
        out_dir=args.out,
        base_seed=args.seed,
        skip_existing=args.skip_existing,
    )


if __name__ == "__main__":
    main()

"""
Drop one or more species from the simulated karyotype and tree files,
then write the reduced versions for reconstruction testing.

Usage:
    python drop_species.py --drop Sp1 Sp2
    python drop_species.py --drop Sp1          # drop one species
"""

import argparse
import re
from pathlib import Path

from ete3 import Tree


def _expand_drop(drop: set, all_names: set) -> set:
    """Also include preWGD variants like Sp2_preWGD for any Sp2 in drop."""
    expanded = set(drop)
    for sp in drop:
        for name in all_names:
            if name == f"{sp}_preWGD":
                expanded.add(name)
    return expanded


def drop_from_karyotypes(src: Path, dst: Path, drop: set):
    """Remove all lines belonging to dropped species (including _preWGD variants)."""
    lines = src.read_text().splitlines()
    # collect all species names present
    all_species = {l[1:].strip() for l in lines if l.startswith(">")}
    drop = _expand_drop(drop, all_species)
    out = []
    skip = False
    for line in lines:
        if line.startswith(">"):
            species = line[1:].strip()
            skip = species in drop
        if not skip:
            out.append(line)
    dst.write_text("\n".join(out) + "\n")
    print(f"Wrote {dst}  (dropped {drop})")


def drop_from_tree(src: Path, dst: Path, drop: set):
    """Prune dropped species from the Newick tree."""
    t = Tree(str(src), format=1)
    for name in drop:
        node = t.search_nodes(name=name)
        if node:
            node[0].detach()
        else:
            print(f"  Warning: {name} not found in tree")
    t.write(outfile=str(dst), format=1)
    print(f"Wrote {dst}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--drop", nargs="+", required=True,
                        help="Species names to drop (e.g. Sp1 Sp2)")
    parser.add_argument("--input-dir", default="output_simulator",
                        help="Source simulator output directory")
    parser.add_argument("--output-dir", default="output_simulator_reduced",
                        help="Destination directory for reduced files")
    args = parser.parse_args()

    drop_base = set(args.drop)
    src = Path(args.input_dir)
    dst = Path(args.output_dir)
    dst.mkdir(parents=True, exist_ok=True)

    # Also drop preWGD counterparts (e.g. Sp2 -> Sp2_preWGD)
    drop = set(drop_base)
    for sp in drop_base:
        drop.add(f"{sp}_preWGD")

    drop_from_karyotypes(
        src / "karyotypes_species_with_outgroup.txt",
        dst / "karyotypes_species_with_outgroup.txt",
        drop,
    )
    drop_from_tree(src / "tree.nwk",               dst / "tree.nwk",               drop - {"Outgroup"})
    drop_from_tree(src / "tree_with_outgroup.nwk",  dst / "tree_with_outgroup.nwk",  drop)

    # Copy files that don't need modification
    for fname in ("karyotypes_true_root.txt",):
        (dst / fname).write_bytes((src / fname).read_bytes())
        print(f"Copied {dst / fname}")

    print(f"\nDone. Now run ancestor_reconstruction.py with:")
    print(f"  input_karyotypes = \"{dst}/karyotypes_species_with_outgroup.txt\"")
    print(f"  input_tree       = \"{dst}/tree.nwk\"")


if __name__ == "__main__":
    main()

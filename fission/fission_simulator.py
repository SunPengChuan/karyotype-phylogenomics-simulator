#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Chromosome Fission Simulator

This script simulates chromosome fission events by breaking chromosomes into two parts.
Unlike misassembly (which cuts a fragment and inserts elsewhere), fission splits one 
chromosome into two independent chromosomes.
"""

import argparse
import csv
import logging
import os
import random
from collections import Counter
from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple


@dataclass
class Chromosome:
    """Represents a chromosome with its genes."""
    species: str
    chr_id: str
    genes: List[str]
    line_index: int

    def label(self) -> str:
        """Return species:chr_id format."""
        return f"{self.species}:{self.chr_id}"


def _setup_logging(log_level: str, log_path: Optional[str]) -> logging.Logger:
    """Setup logging with both console and optional file output."""
    logger = logging.getLogger("fission_simulator")
    logger.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    logger.handlers.clear()

    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    if log_path:
        file_handler = logging.FileHandler(log_path, encoding="utf-8")
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def _parse_karyotypes(path: str) -> List[Chromosome]:
    """Parse karyotypes file into Chromosome objects."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Input file not found: {path}")

    chromosomes: List[Chromosome] = []
    with open(path, "r", encoding="utf-8") as f:
        lines = [line.rstrip("\n") for line in f]

    if len(lines) < 3:
        raise ValueError("Input karyotypes file is too short.")

    for idx, raw in enumerate(lines[2:], start=3):
        line = raw.strip()
        if not line:
            continue
        if set(line) == {"-"}:
            continue
        parts = line.split()
        if len(parts) < 3:
            raise ValueError(f"Invalid karyotypes line {idx}: expected >=3 columns, got {len(parts)}")
        species, chr_id, *genes = parts
        if not genes:
            raise ValueError(f"Invalid karyotypes line {idx}: empty gene list")
        chromosomes.append(Chromosome(species=species, chr_id=chr_id, genes=genes, line_index=idx))

    if not chromosomes:
        raise ValueError("No chromosome records parsed from input.")
    return chromosomes


def _write_karyotypes(path: str, chromosomes: Sequence[Chromosome]) -> None:
    """Write chromosomes to karyotypes file."""
    header1 = "Species         Chr             Genes"
    header2 = "-" * 80
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(header1 + "\n")
        f.write(header2 + "\n")
        for chrom in chromosomes:
            genes_str = " ".join(chrom.genes)
            f.write(f"{chrom.species:<12}{chrom.chr_id:<16}{genes_str}\n")


def _generate_split_id(original_id: str, suffix: str) -> str:
    """
    Generate chromosome ID after fission.
    
    Examples:
        "1" -> "1a", "1b" or "1(1)", "1(2)"
        "A" -> "Aa", "Ab" or "A(1)", "A(2)"
    """
    return f"{original_id}{suffix}"


def _pick_fission_point(
    rng: random.Random,
    genes: Sequence[str],
    min_block_len: int,
) -> int:
    """
    Pick a random position to split chromosome into two.
    
    Args:
        rng: Random number generator
        genes: Gene sequence
        min_block_len: Minimum genes required in each resulting chromosome
    
    Returns:
        Split position (0-based, split before this position)
    
    Example:
        genes = [A, B, C, D, E], min_block_len=2
        valid split positions: 2 or 3
        if split_pos=2: left=[A,B], right=[C,D,E]
    """
    n = len(genes)
    if n < 2 * min_block_len:
        raise ValueError(f"Chromosome too short to split (needs >= {2 * min_block_len} genes, has {n})")
    
    # Valid split positions: from min_block_len to n - min_block_len
    valid_starts = list(range(min_block_len, n - min_block_len + 1))
    if not valid_starts:
        raise ValueError(f"No valid split positions with min_block_len={min_block_len}")
    
    return rng.choice(valid_starts)


def simulate_fission(
    chromosomes: List[Chromosome],
    seed: int,
    fission_rate: float,
    max_fissions: Optional[int],
    min_block_len: int,
    id_suffix_mode: str,
    logger: logging.Logger,
) -> Tuple[List[dict], int]:
    """
    Simulate chromosome fission events.
    
    Args:
        chromosomes: List of chromosomes to process
        seed: Random seed for reproducibility
        fission_rate: Fraction of chromosomes to split (0-1)
        max_fissions: Optional cap on number of fissions
        min_block_len: Minimum genes in each resulting chromosome
        id_suffix_mode: How to name split chromosomes ('letter' or 'number')
        logger: Logger instance
    
    Returns:
        Tuple of (report_rows, number_of_fissions_applied)
    """
    rng = random.Random(seed)

    # Only split species chromosomes (Sp*) with sufficient length
    eligible_indices = [
        i for i, c in enumerate(chromosomes) 
        if len(c.genes) >= 2 * min_block_len and c.species.startswith("Sp")
    ]
    
    if not eligible_indices:
        logger.info("No eligible chromosomes for fission (need Sp* with >=%d genes).", 2 * min_block_len)
        return [], 0

    n_eligible = len(eligible_indices)
    k = int(n_eligible * fission_rate)
    if max_fissions is not None and max_fissions >= 0:
        k = min(k, max_fissions)
    k = min(k, n_eligible)
    
    if k <= 0:
        logger.info("No fission events to apply (rate=%.3f, max=%s).", fission_rate, max_fissions)
        return [], 0

    # Randomly select chromosomes to split
    selected_indices = rng.sample(eligible_indices, k)
    
    report_rows: List[dict] = []
    new_chromosomes: List[Chromosome] = []
    indices_to_remove: set = set()
    
    logger.info(
        "Total chromosomes=%d; eligible=%d; planned fissions=%d; seed=%d",
        len(chromosomes), n_eligible, k, seed
    )

    for idx in selected_indices:
        chrom = chromosomes[idx]
        original_len = len(chrom.genes)
        
        try:
            split_pos = _pick_fission_point(rng, chrom.genes, min_block_len)
        except ValueError as e:
            logger.warning("Skipping %s: %s", chrom.label(), str(e))
            continue

        # Split chromosome
        left_genes = chrom.genes[:split_pos]
        right_genes = chrom.genes[split_pos:]
        
        # Generate new chromosome IDs
        if id_suffix_mode == "letter":
            left_suffix = "a"
            right_suffix = "b"
        else:  # number mode
            left_suffix = "(1)"
            right_suffix = "(2)"
        
        left_id = _generate_split_id(chrom.chr_id, left_suffix)
        right_id = _generate_split_id(chrom.chr_id, right_suffix)
        
        # Create new chromosomes
        left_chrom = Chromosome(
            species=chrom.species,
            chr_id=left_id,
            genes=left_genes,
            line_index=chrom.line_index
        )
        right_chrom = Chromosome(
            species=chrom.species,
            chr_id=right_id,
            genes=right_genes,
            line_index=chrom.line_index
        )
        
        # Mark original for removal and add new ones
        indices_to_remove.add(idx)
        new_chromosomes.extend([left_chrom, right_chrom])
        
        # Record the event
        row = dict(
            seed=seed,
            species=chrom.species,
            original_chr=chrom.chr_id,
            original_len=original_len,
            split_pos_0based=split_pos,
            left_chr=left_id,
            left_len=len(left_genes),
            right_chr=right_id,
            right_len=len(right_genes),
        )
        report_rows.append(row)
        
        logger.info(
            "Fission: %s (%d genes) -> %s (%d) + %s (%d) @ position %d",
            chrom.label(),
            original_len,
            left_id,
            len(left_genes),
            right_id,
            len(right_genes),
            split_pos
        )

    # Build final chromosome list: remove split ones and add new ones
    final_chromosomes = [c for i, c in enumerate(chromosomes) if i not in indices_to_remove]
    final_chromosomes.extend(new_chromosomes)
    
    # Update the original list in-place
    chromosomes.clear()
    chromosomes.extend(final_chromosomes)

    return report_rows, len(report_rows)


def _write_report(path: str, rows: Sequence[dict]) -> None:
    """Write fission report to TSV file."""
    fieldnames = [
        "Species",
        "OriginalChr",
        "OriginalLen",
        "SplitPos",
        "LeftChr",
        "LeftLen",
        "RightChr",
        "RightLen",
    ]
    
    simplified_rows = []
    for r in rows:
        simplified_rows.append({
            "Species": r["species"],
            "OriginalChr": r["original_chr"],
            "OriginalLen": r["original_len"],
            "SplitPos": r["split_pos_0based"],
            "LeftChr": r["left_chr"],
            "LeftLen": r["left_len"],
            "RightChr": r["right_chr"],
            "RightLen": r["right_len"],
        })
    
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for row in simplified_rows:
            w.writerow(row)


def _read_report(path: str) -> List[dict]:
    """Read fission report from TSV file."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Report file not found: {path}")
    with open(path, "r", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _verify_rows(rows: Sequence[dict]) -> dict:
    """Verify and summarize fission report."""
    if not rows:
        return dict(rows=0, unique_species=0, total_genes_before=0, total_genes_after=0)
    
    species_counter = Counter(r["species"] for r in rows)
    total_before = sum(int(r["original_len"]) for r in rows)
    total_left = sum(int(r["left_len"]) for r in rows)
    total_right = sum(int(r["right_len"]) for r in rows)
    
    return dict(
        rows=len(rows),
        unique_species=len(species_counter),
        total_genes_before=total_before,
        total_genes_after=total_left + total_right,
        genes_preserved=(total_before == total_left + total_right),
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Simulate chromosome fission by breaking chromosomes into two parts."
    )
    parser.add_argument(
        "--input",
        default="output_simulator/karyotypes_species_with_outgroup.txt",
        help="Input karyotypes file path (default: output_simulator/karyotypes_species_with_outgroup.txt)"
    )
    parser.add_argument(
        "--output-karyotypes",
        required=True,
        help="Output karyotypes file path"
    )
    parser.add_argument(
        "--report",
        required=True,
        help="Output fission report TSV path"
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed (default: 42)"
    )
    parser.add_argument(
        "--fission-rate",
        type=float,
        default=0.01,
        help="Fraction of chromosomes to split (default: 0.01, i.e., 10%%)"
    )
    parser.add_argument(
        "--max-fissions",
        type=int,
        default=None,
        help="Optional cap on number of fission events (default: no limit)"
    )
    parser.add_argument(
        "--min-block-len",
        type=int,
        default=20,
        help="Minimum genes in each resulting chromosome (default: 20)"
    )
    parser.add_argument(
        "--id-suffix-mode",
        choices=["letter", "number"],
        default="letter",
        help="How to name split chromosomes: 'letter' (1->1a,1b) or 'number' (1->1(1),1(2))"
    )
    parser.add_argument(
        "--log",
        default=None,
        help="Optional log file path"
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        help="Log level: DEBUG/INFO/WARNING/ERROR (default: INFO)"
    )
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="Only verify an existing report TSV and exit"
    )

    args = parser.parse_args(argv)
    logger = _setup_logging(args.log_level, args.log)

    try:
        if args.verify_only:
            rows = _read_report(args.report)
            summary = _verify_rows(rows)
            logger.info("Verify summary: %s", summary)
            return 0

        # Parse input karyotypes
        chromosomes = _parse_karyotypes(args.input)
        logger.info("Loaded %d chromosomes from %s", len(chromosomes), args.input)
        
        # Simulate fission
        report_rows, applied = simulate_fission(
            chromosomes=chromosomes,
            seed=args.seed,
            fission_rate=args.fission_rate,
            max_fissions=args.max_fissions,
            min_block_len=args.min_block_len,
            id_suffix_mode=args.id_suffix_mode,
            logger=logger,
        )
        
        # Write outputs
        _write_karyotypes(args.output_karyotypes, chromosomes)
        _write_report(args.report, report_rows)
        
        # Verify and summarize
        summary = _verify_rows(report_rows)
        logger.info("Verify summary: %s", summary)
        logger.info(
            "Done. Applied fissions=%d; output=%s; report=%s",
            applied, args.output_karyotypes, args.report
        )
        
        return 0
        
    except Exception as e:
        logger.exception("Failed: %s", str(e))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

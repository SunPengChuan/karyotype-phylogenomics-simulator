import argparse
import collections
import csv
import logging
import os
import random
from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence, Tuple


@dataclass
class Chromosome:
    species: str
    chr_id: str
    genes: List[str]
    line_index: int

    def label(self) -> str:
        return f"{self.species}:{self.chr_id}"


def _setup_logging(log_level: str, log_path: Optional[str]) -> logging.Logger:
    logger = logging.getLogger("misassembly_simulator")
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
    header1 = "Species         Chr             Genes"
    header2 = "-" * 80
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(header1 + "\n")
        f.write(header2 + "\n")
        for chrom in chromosomes:
            genes_str = " ".join(chrom.genes)
            f.write(f"{chrom.species:<12}{chrom.chr_id:<16}{genes_str}\n")


def _pick_fragment(
    rng: random.Random,
    genes: Sequence[str],
    max_fragment_ratio: float,
) -> Tuple[int, int]:
    if not (0.0 < max_fragment_ratio <= 0.5):
        raise ValueError("--max-fragment-ratio must be within (0, 0.5]")
    l_total = len(genes)
    if l_total < 2:
        raise ValueError("Source chromosome too short to cut (needs at least 2 genes).")

    max_len = max(1, int(l_total * max_fragment_ratio))
    max_len = min(max_len, l_total - 1)
    frag_len = rng.randint(1, max_len)
    start = rng.randint(0, l_total - frag_len)
    end = start + frag_len - 1
    return start, end


def simulate_misassembly(
    chromosomes: List[Chromosome],
    seed: int,
    error_rate: float,
    max_errors: Optional[int],
    max_fragment_ratio: float,
    logger: logging.Logger,
) -> Tuple[List[dict], int]:
    rng = random.Random(seed)

    eligible_indices = [i for i, c in enumerate(chromosomes) if len(c.genes) >= 2 and c.species.startswith("Sp")]
    if len(eligible_indices) < 2:
        logger.info("Not enough eligible chromosomes (need >=2 Sp chromosomes with >=2 genes). No misassemblies applied.")
        return [], 0

    n_eligible = len(eligible_indices)
    k = int(n_eligible * error_rate)
    if max_errors is not None and max_errors >= 0:
        k = min(k, max_errors)
    k = min(k, n_eligible // 2)
    if k <= 0:
        return [], 0

    sources = rng.sample(eligible_indices, k)
    remaining = [i for i in eligible_indices if i not in set(sources)]
    targets = rng.sample(remaining, k) if len(remaining) >= k else []

    report_rows: List[dict] = []
    logger.info("Total chromosomes=%d; eligible=%d; planned misassemblies=%d; seed=%d", len(chromosomes), n_eligible, len(sources), seed)

    for src_i, tgt_i in zip(sources, targets):
        src = chromosomes[src_i]
        tgt = chromosomes[tgt_i]
        if src_i == tgt_i:
            continue

        src_len_before = len(src.genes)
        tgt_len_before = len(tgt.genes)
        try:
            start, end = _pick_fragment(rng, src.genes, max_fragment_ratio)
        except ValueError as e:
            logger.warning("Skipping source %s: %s", src.label(), str(e))
            continue

        fragment = src.genes[start : end + 1]
        del src.genes[start : end + 1]

        insert_pos = rng.randint(0, len(tgt.genes))
        tgt.genes[insert_pos:insert_pos] = fragment

        row = dict(
            seed=seed,
            source_species=src.species,
            source_chr=src.chr_id,
            source_len_before=src_len_before,
            cut_start_0based=start,
            cut_end_0based=end,
            cut_len=len(fragment),
            target_species=tgt.species,
            target_chr=tgt.chr_id,
            target_len_before=tgt_len_before,
            insert_pos_0based=insert_pos,
        )
        report_rows.append(row)
        logger.info(
            "Misassembly: %s cut[%d,%d] len=%d -> %s insert@%d",
            src.label(),
            start,
            end,
            len(fragment),
            tgt.label(),
            insert_pos,
        )

    return report_rows, len(report_rows)


def _write_report(path: str, rows: Sequence[dict]) -> None:
    simplified_rows = []
    for r in rows:
        simplified_rows.append({
            "Source": f"{r['source_species']}:{r['source_chr']}",
            "Cut": f"{r['cut_start_0based']}-{r['cut_end_0based']}",
            "Length": r['cut_len'],
            "Target": f"{r['target_species']}:{r['target_chr']}",
            "Insert": r['insert_pos_0based'],
        })
    fieldnames = ["Source", "Cut", "Length", "Target", "Insert"]
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for row in simplified_rows:
            w.writerow(row)


def _read_report(path: str) -> List[dict]:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Report file not found: {path}")
    with open(path, "r", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _verify_rows(rows: Sequence[dict]) -> dict:
    sources = [(r.get("source_species"), r.get("source_chr")) for r in rows]
    targets = [(r.get("target_species"), r.get("target_chr")) for r in rows]
    dup_sources = sum(1 for _, c in collections.Counter(sources).items() if c > 1)
    dup_targets = sum(1 for _, c in collections.Counter(targets).items() if c > 1)
    return dict(
        rows=len(rows),
        unique_sources=len(set(sources)),
        unique_targets=len(set(targets)),
        dup_sources=dup_sources,
        dup_targets=dup_targets,
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Simulate chromosome misassembly errors by cutting fragments and inserting into other chromosomes.")
    parser.add_argument(
        "--input",
        default=r"d:\git\cell\respone1\ceshi\output_simulator\karyotypes_species_with_outgroup.txt",
        help="Input karyotypes_species.txt path",
    )
    parser.add_argument("--output-karyotypes", required=True, help="Output karyotypes file path")
    parser.add_argument("--report", required=True, help="Output misassembly report TSV path")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--error-rate", type=float, default=0.01, help="Fraction of chromosomes to misassemble (<=0.1)")
    parser.add_argument("--max-errors", type=int, default=None, help="Optional absolute cap on misassemblies (0 disables)")
    parser.add_argument("--max-fragment-ratio", type=float, default=0.5, help="Max fragment length ratio (<=0.5)")
    parser.add_argument("--log", default=None, help="Optional log file path")
    parser.add_argument("--log-level", default="INFO", help="Log level: DEBUG/INFO/WARNING/ERROR")
    parser.add_argument("--verify-only", action="store_true", help="Only verify an existing report TSV and exit")

    args = parser.parse_args(argv)
    logger = _setup_logging(args.log_level, args.log)

    try:
        if args.verify_only:
            rows = _read_report(args.report)
            summary = _verify_rows(rows)
            logger.info("Verify summary: %s", summary)
            return 0

        chromosomes = _parse_karyotypes(args.input)
        report_rows, applied = simulate_misassembly(
            chromosomes=chromosomes,
            seed=args.seed,
            error_rate=args.error_rate,
            max_errors=args.max_errors,
            max_fragment_ratio=args.max_fragment_ratio,
            logger=logger,
        )
        _write_karyotypes(args.output_karyotypes, chromosomes)
        _write_report(args.report, report_rows)
        summary = _verify_rows(report_rows)
        logger.info("Verify summary: %s", summary)
        logger.info("Done. Applied misassemblies=%d; output=%s; report=%s", applied, args.output_karyotypes, args.report)
        return 0
    except Exception as e:
        logger.exception("Failed: %s", str(e))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

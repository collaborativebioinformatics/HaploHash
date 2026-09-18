#!/usr/bin/env python3
"""Aggregate AlphaGenome Atlas AVI scores into hg38 haploblocks.

The static AVI download contains one unsigned prioritisation score per alternate
allele. This script streams all rows overlapping the supplied, non-overlapping
BED intervals and reports block-level summaries without loading a chromosome
into memory.

AVI PHRED is a genome-wide rank, not a signed molecular effect. Consequently,
the summaries retain extremes and PHRED >= 20 hit counts rather than treating
the sum of PHRED values as an additive biological effect.

Examples:
  python3 score_blocks.py blocks.bed avi.tsv.gz --inspect
  python3 score_blocks.py blocks.bed avi.tsv.gz block_scores.tsv
  python3 score_blocks.py blocks.bed avi.tsv.gz block_scores.tsv \
      --phred-column 6

Column numbers supplied on the command line are 1-based. A matching .tbi index
must sit beside a local AVI file or be available beside a range-queryable URL.
"""

from __future__ import annotations

import argparse
import heapq
import math
import re
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator, Sequence


PHRED_THRESHOLD = 20.0  # Atlas paper: top 1% genome-wide
DEFAULT_TOP_K = 10


@dataclass(frozen=True)
class Block:
    chrom: str
    start: int
    end: int
    block_id: str

    @property
    def length_bp(self) -> int:
        return self.end - self.start


@dataclass
class BlockSummary:
    top_k: int
    n_positions: int = 0
    n_variant_alleles: int = 0
    n_scored: int = 0
    n_missing_avi: int = 0
    n_top1: int = 0
    sum_top1: float = 0.0
    sum_all: float = 0.0
    maximum: float | None = None
    _last_position: int | None = None
    _top_scores: list[float] = field(default_factory=list)

    def add(self, position: int, phred: float | None) -> None:
        """Add one alternate-allele row at a 1-based genomic position."""
        self.n_variant_alleles += 1
        if position != self._last_position:
            self.n_positions += 1
            self._last_position = position

        if phred is None or not math.isfinite(phred):
            self.n_missing_avi += 1
            return

        self.n_scored += 1
        self.sum_all += phred
        self.maximum = phred if self.maximum is None else max(self.maximum, phred)
        if len(self._top_scores) < self.top_k:
            heapq.heappush(self._top_scores, phred)
        elif phred > self._top_scores[0]:
            heapq.heapreplace(self._top_scores, phred)

        if phred >= PHRED_THRESHOLD:
            self.n_top1 += 1
            self.sum_top1 += phred

    def values(self, length_bp: int) -> list[str]:
        top_k_mean = (
            sum(self._top_scores) / len(self._top_scores)
            if self._top_scores
            else None
        )
        top1_fraction = self.n_top1 / self.n_scored if self.n_scored else None
        top1_per_kb = (
            self.n_top1 / (length_bp / 1000.0) if length_bp > 0 else None
        )
        top1_mean = self.sum_top1 / self.n_top1 if self.n_top1 else None
        avi_mean = self.sum_all / self.n_scored if self.n_scored else None
        return [
            str(self.n_positions),
            str(self.n_variant_alleles),
            str(self.n_scored),
            str(self.n_missing_avi),
            format_number(self.maximum),
            format_number(top_k_mean),
            str(self.n_top1),
            format_number(top1_fraction),
            format_number(top1_per_kb),
            format_number(top1_mean),
            format_number(avi_mean),
        ]


def format_number(value: float | None) -> str:
    return "" if value is None else f"{value:.8g}"


def canonical_chrom(chrom: str) -> str:
    chrom = chrom.strip()
    return chrom if chrom.lower().startswith("chr") else f"chr{chrom}"


def normalise_column_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", name.lstrip("#").lower())


def parse_column_spec(
    spec: str | None,
    header: Sequence[str] | None,
    aliases: Sequence[str],
    label: str,
    fallback: int | None = None,
) -> int:
    """Resolve a 1-based integer or header name to a 0-based index."""
    if spec is not None:
        try:
            index = int(spec) - 1
        except ValueError:
            if header is None:
                raise ValueError(
                    f"--{label}-column={spec!r} is a name, but no header was found"
                ) from None
            wanted = normalise_column_name(spec)
            matches = [
                i
                for i, column in enumerate(header)
                if normalise_column_name(column) == wanted
            ]
            if len(matches) != 1:
                raise ValueError(
                    f"could not resolve --{label}-column={spec!r} in header {list(header)!r}"
                )
            return matches[0]
        if index < 0:
            raise ValueError(f"--{label}-column must be a positive 1-based number")
        if header is not None and index >= len(header):
            raise ValueError(
                f"--{label}-column={spec} exceeds the {len(header)} header columns"
            )
        return index

    if header is not None:
        normalised = [normalise_column_name(column) for column in header]
        for alias in aliases:
            alias_normalised = normalise_column_name(alias)
            if alias_normalised in normalised:
                return normalised.index(alias_normalised)

    if fallback is not None:
        return fallback
    raise ValueError(
        f"could not infer the {label} column; pass --{label}-column NAME_OR_NUMBER"
    )


def run_tabix(
    tabix: str, arguments: Sequence[str], *, capture_output: bool = True
) -> subprocess.CompletedProcess[str]:
    try:
        return subprocess.run(
            [tabix, *arguments],
            check=True,
            capture_output=capture_output,
            text=True,
        )
    except FileNotFoundError:
        raise RuntimeError(
            f"tabix executable not found: {tabix!r}; install htslib/tabix or pass --tabix"
        ) from None
    except subprocess.CalledProcessError as error:
        detail = (error.stderr or error.stdout or "").strip()
        raise RuntimeError(f"tabix failed: {detail or error}") from error


def read_tabix_header(tabix: str, avi_path: str) -> list[str] | None:
    result = run_tabix(tabix, ["-H", avi_path])
    header_lines = [line for line in result.stdout.splitlines() if line.startswith("#")]
    if not header_lines:
        return None
    return header_lines[-1].lstrip("#").split("\t")


def read_tabix_contigs(tabix: str, avi_path: str) -> list[str]:
    result = run_tabix(tabix, ["-l", avi_path])
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def detect_contig_style(contigs: Sequence[str]) -> str:
    if not contigs:
        raise ValueError("the Tabix index contains no contigs")
    has_chr = [contig.lower().startswith("chr") for contig in contigs]
    if all(has_chr):
        return "chr"
    if not any(has_chr):
        return "bare"
    raise ValueError(
        "the Tabix index mixes chr-prefixed and bare contigs; pass --contig-style"
    )


def tabix_chrom(chrom: str, style: str) -> str:
    canonical = canonical_chrom(chrom)
    return canonical if style == "chr" else canonical[3:]


def load_blocks(path: Path) -> list[Block]:
    blocks: list[Block] = []
    with path.open() as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError(
                    f"{path}:{line_number}: expected at least 4 tab-separated columns"
                )
            chrom, start_text, end_text, block_id = fields[:4]
            try:
                start, end = int(start_text), int(end_text)
            except ValueError:
                raise ValueError(
                    f"{path}:{line_number}: start/end must be integers"
                ) from None
            if start < 0 or end <= start:
                raise ValueError(
                    f"{path}:{line_number}: invalid half-open interval {start}-{end}"
                )
            blocks.append(Block(canonical_chrom(chrom), start, end, block_id))

    if not blocks:
        raise ValueError(f"no blocks found in {path}")

    seen_ids: set[str] = set()
    last_by_chrom: dict[str, Block] = {}
    for block in blocks:
        if block.block_id in seen_ids:
            raise ValueError(f"duplicate block id: {block.block_id}")
        seen_ids.add(block.block_id)
        previous = last_by_chrom.get(block.chrom)
        if previous is not None:
            if block.start < previous.start:
                raise ValueError(f"blocks are not sorted within {block.chrom}")
            if block.start < previous.end:
                raise ValueError(
                    f"overlapping blocks: {previous.block_id} and {block.block_id}"
                )
        last_by_chrom[block.chrom] = block
    return blocks


def write_tabix_regions(blocks: Sequence[Block], style: str) -> Path:
    handle = tempfile.NamedTemporaryFile(
        mode="w", prefix="haploblocks-", suffix=".bed", delete=False
    )
    path = Path(handle.name)
    with handle:
        for block in blocks:
            handle.write(
                f"{tabix_chrom(block.chrom, style)}\t{block.start}\t{block.end}\n"
            )
    return path


def stream_tabix_rows(
    tabix: str, avi_path: str, regions_path: Path
) -> Iterator[str]:
    try:
        process = subprocess.Popen(
            [tabix, "-R", str(regions_path), avi_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except FileNotFoundError:
        raise RuntimeError(f"tabix executable not found: {tabix!r}") from None

    assert process.stdout is not None
    for line in process.stdout:
        if line and not line.startswith("#"):
            yield line.rstrip("\n")
    stderr = process.stderr.read() if process.stderr is not None else ""
    return_code = process.wait()
    if return_code:
        raise RuntimeError(f"tabix failed ({return_code}): {stderr.strip()}")


def block_lookup(blocks: Sequence[Block]) -> tuple[dict[str, list[Block]], dict[str, int]]:
    by_chrom: dict[str, list[Block]] = {}
    for block in blocks:
        by_chrom.setdefault(block.chrom, []).append(block)
    return by_chrom, {chrom: 0 for chrom in by_chrom}


def aggregate_rows(
    rows: Iterable[str],
    blocks: Sequence[Block],
    *,
    chrom_column: int,
    position_column: int,
    phred_column: int,
    top_k: int,
) -> tuple[list[BlockSummary], int]:
    """Aggregate sorted Tabix rows and return summaries plus unmatched rows."""
    summaries = [BlockSummary(top_k=top_k) for _ in blocks]
    summary_by_id = {
        block.block_id: summary for block, summary in zip(blocks, summaries, strict=True)
    }
    by_chrom, cursors = block_lookup(blocks)
    unmatched = 0

    required_column = max(chrom_column, position_column, phred_column)
    for line_number, line in enumerate(rows, 1):
        fields = line.split("\t")
        if len(fields) <= required_column:
            raise ValueError(
                f"Tabix row {line_number} has {len(fields)} columns; need column "
                f"{required_column + 1}: {line[:200]!r}"
            )
        chrom = canonical_chrom(fields[chrom_column])
        try:
            position = int(fields[position_column])
        except ValueError:
            raise ValueError(
                f"Tabix row {line_number} has a non-integer position: "
                f"{fields[position_column]!r}"
            ) from None

        chrom_blocks = by_chrom.get(chrom)
        if not chrom_blocks:
            unmatched += 1
            continue
        cursor = cursors[chrom]
        position_zero_based = position - 1
        while (
            cursor < len(chrom_blocks)
            and chrom_blocks[cursor].end <= position_zero_based
        ):
            cursor += 1
        cursors[chrom] = cursor
        if cursor >= len(chrom_blocks):
            unmatched += 1
            continue
        block = chrom_blocks[cursor]
        if not (block.start <= position_zero_based < block.end):
            unmatched += 1
            continue

        phred_text = fields[phred_column].strip()
        try:
            phred = (
                float(phred_text)
                if phred_text not in {"", ".", "NA", "NaN"}
                else None
            )
        except ValueError:
            phred = None
        summary_by_id[block.block_id].add(position, phred)

    return summaries, unmatched


def inspect_file(
    args: argparse.Namespace,
    blocks: Sequence[Block],
    header: Sequence[str] | None,
    contig_style: str,
) -> int:
    first = blocks[0]
    region = (
        f"{tabix_chrom(first.chrom, contig_style)}:"
        f"{first.start + 1}-{first.end}"
    )
    result = run_tabix(args.tabix, [args.avi_path, region])
    rows = [line for line in result.stdout.splitlines() if line][:5]
    print(f"header: {list(header) if header else '(none)'}")
    print(f"first block: {first.block_id} ({region}); showing {len(rows)} rows")
    for row in rows:
        fields = row.split("\t")
        print(f"  {len(fields)} columns: {fields}")

    try:
        indices = resolve_columns(args, header)
    except ValueError as error:
        print(f"\ncolumn resolution: {error}")
        print("Pass explicit 1-based column numbers before running the full aggregation.")
        return 0
    print(
        "\nresolved columns (1-based): "
        f"chrom={indices[0] + 1}, position={indices[1] + 1}, "
        f"phred={indices[2] + 1}"
    )
    return 0


def resolve_columns(
    args: argparse.Namespace, header: Sequence[str] | None
) -> tuple[int, int, int]:
    chrom_column = parse_column_spec(
        args.chrom_column,
        header,
        aliases=("chrom", "chr", "chromosome", "seqname"),
        label="chrom",
        fallback=0,
    )
    position_column = parse_column_spec(
        args.position_column,
        header,
        aliases=("pos", "position", "bp"),
        label="position",
        fallback=1,
    )
    phred_column = parse_column_spec(
        args.phred_column,
        header,
        aliases=("phred", "avi_phred", "phred_score", "avi_phred_score"),
        label="phred",
    )
    return chrom_column, position_column, phred_column


def write_output(
    output_path: Path,
    blocks: Sequence[Block],
    summaries: Sequence[BlockSummary],
    top_k: int,
) -> None:
    with output_path.open("w") as output:
        output.write(
            "block_id\tchrom\tstart\tend\tlength_bp\tn_positions\t"
            "n_variant_alleles\tn_scored\tn_missing_avi\tavi_max\t"
            f"avi_top{top_k}_mean\tavi_top1_count\tavi_top1_fraction\t"
            "avi_top1_per_kb\tavi_top1_mean\tavi_mean\n"
        )
        for block, summary in zip(blocks, summaries, strict=True):
            fixed = [
                block.block_id,
                block.chrom,
                str(block.start),
                str(block.end),
                str(block.length_bp),
            ]
            output.write("\t".join([*fixed, *summary.values(block.length_bp)]) + "\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("blocks_bed", type=Path)
    parser.add_argument("avi_path", help="local path or HTTP(S) URL to AVI .tsv.gz")
    parser.add_argument(
        "output_tsv", nargs="?", type=Path, default=Path("block_scores.tsv")
    )
    parser.add_argument(
        "--inspect",
        action="store_true",
        help="show the header, sample rows and resolved columns, then exit",
    )
    parser.add_argument("--tabix", default="tabix", help="tabix executable")
    parser.add_argument(
        "--chrom-column", help="chromosome column name or 1-based number"
    )
    parser.add_argument(
        "--position-column", help="position column name or 1-based number"
    )
    parser.add_argument(
        "--phred-column", help="AVI PHRED column name or 1-based number"
    )
    parser.add_argument(
        "--contig-style",
        choices=("auto", "chr", "bare"),
        default="auto",
        help="contig naming in the AVI index (default: detect with tabix -l)",
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=DEFAULT_TOP_K,
        help=f"number of strongest alleles in the top-k mean (default: {DEFAULT_TOP_K})",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.top_k <= 0:
        raise ValueError("--top-k must be positive")

    blocks = load_blocks(args.blocks_bed)
    header = read_tabix_header(args.tabix, args.avi_path)
    contig_style = args.contig_style
    if contig_style == "auto":
        contig_style = detect_contig_style(
            read_tabix_contigs(args.tabix, args.avi_path)
        )

    if args.inspect:
        return inspect_file(args, blocks, header, contig_style)

    chrom_column, position_column, phred_column = resolve_columns(args, header)
    regions_path = write_tabix_regions(blocks, contig_style)
    started = time.time()
    try:
        summaries, unmatched = aggregate_rows(
            stream_tabix_rows(args.tabix, args.avi_path, regions_path),
            blocks,
            chrom_column=chrom_column,
            position_column=position_column,
            phred_column=phred_column,
            top_k=args.top_k,
        )
    finally:
        regions_path.unlink(missing_ok=True)

    write_output(args.output_tsv, blocks, summaries, args.top_k)
    elapsed = time.time() - started
    total_rows = sum(summary.n_variant_alleles for summary in summaries)
    print(
        f"wrote {args.output_tsv}: {len(blocks)} blocks, "
        f"{total_rows} variant alleles in {elapsed:.1f}s",
        file=sys.stderr,
    )
    if unmatched:
        print(
            f"WARNING: skipped {unmatched} Tabix rows outside the supplied blocks",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (RuntimeError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2) from error

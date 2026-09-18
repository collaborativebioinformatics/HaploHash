#!/usr/bin/env python3
"""Filter the chr22 pilot to observed SNVs, calculate AF/MAF, and join AVI blocks."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import Counter
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from build_variant_parquet import INDEL_COLUMNS, SNV_COLUMNS, pilot_blocks
from score_blocks import Block, canonical_chrom


HERE = Path(__file__).resolve().parent
BLOCK_INTS = (
    "start", "end", "length_bp", "n_positions", "n_variant_alleles",
    "n_scored", "n_missing_avi", "avi_top1_count",
)
ADDED_FIELDS = [
    pa.field("block_key", pa.string()),
    pa.field("block_id_avi", pa.string()),
    pa.field("AF_from_counts", pa.float64()),
    pa.field("MAF_from_counts", pa.float64()),
    pa.field("has_avi_score", pa.bool_()),
]


def coordinate_key(chrom: str, start: int, end: int) -> tuple[str, int, int]:
    return canonical_chrom(chrom), start, end


def integer(value: object, label: str) -> int:
    try:
        result = int(value)
        if isinstance(value, bool) or float(value) != result:
            raise ValueError
        return result
    except (TypeError, ValueError, OverflowError):
        raise ValueError(f"{label}: expected an integer, got {value!r}") from None


def require_columns(columns, required, label: str) -> None:
    missing = sorted(set(required) - set(columns))
    if missing:
        raise ValueError(f"{label}: missing columns {missing}")


def load_block_scores(path: Path) -> list[dict]:
    with path.open() as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        require_columns(reader.fieldnames or [], ["block_id", "chrom", *BLOCK_INTS], str(path))
        rows = list(reader)
    for row in rows:
        for column in BLOCK_INTS:
            row[column] = integer(row[column], f"{row['block_id']} {column}")
    return rows


def select_blocks(rows: list[dict], expected: list[Block]) -> dict[tuple, dict]:
    """Require exactly one AVI row for every pilot interval, including empty ones."""
    index = {}
    for row in rows:
        key = coordinate_key(row["chrom"], row["start"], row["end"])
        if key in index:
            raise ValueError(f"duplicate block coordinates: {key}")
        index[key] = row
    selected = {}
    for block in expected:
        key = coordinate_key(block.chrom, block.start, block.end)
        if key not in index:
            raise ValueError(f"missing AVI block: {key}")
        row = index[key]
        if row["length_bp"] != block.end - block.start:
            raise ValueError(f"block length disagrees with coordinates: {key}")
        if not (
            0 <= row["n_positions"] <= row["length_bp"]
            and 0 <= row["avi_top1_count"] <= row["n_scored"] <= row["n_variant_alleles"]
            and row["n_missing_avi"] >= 0
            and row["n_scored"] + row["n_missing_avi"] == row["n_variant_alleles"]
            and row["n_variant_alleles"] <= 3 * row["n_positions"]
        ):
            raise ValueError(f"inconsistent AVI allele counts: {key}")
        selected[key] = row
    return selected


def validate_variant(source: dict, blocks: dict, *, snv: bool) -> dict:
    row = dict(source)
    label = str(row.get("variant_id", "variant"))
    require_columns(row, SNV_COLUMNS if snv else INDEL_COLUMNS, label)
    for column in ("pos", "AC", "AN", "block_start", "block_end"):
        row[column] = integer(row[column], f"{label} {column}")
    if not isinstance(row["chrom"], str) or not row["chrom"].strip():
        raise ValueError(f"{label}: missing chromosome")
    row["chrom"] = canonical_chrom(row["chrom"])
    key = coordinate_key(row["chrom"], row["block_start"], row["block_end"])
    if key not in blocks:
        raise ValueError(f"{label}: assigned block is not a pilot interval: {key}")
    if not row["block_start"] <= row["pos"] - 1 < row["block_end"]:
        raise ValueError(f"{label}: position outside assigned half-open block")
    if not (row["AN"] > 0 and 0 <= row["AC"] <= row["AN"]):
        raise ValueError(f"{label}: invalid AC/AN: {row['AC']}/{row['AN']}")
    try:
        row["AF"] = float(row["AF"])
    except (TypeError, ValueError):
        raise ValueError(f"{label}: invalid reported AF") from None
    if not math.isfinite(row["AF"]) or not 0 <= row["AF"] <= 1:
        raise ValueError(f"{label}: reported AF outside [0, 1]")
    ref, alt = row["ref"], row["alt"]
    if any(not isinstance(a, str) or not a or set(a) - set("ACGT") for a in (ref, alt)):
        raise ValueError(f"{label}: REF/ALT must contain only uppercase A/C/G/T")
    if ref == alt or (snv and (len(ref) != 1 or len(alt) != 1)):
        raise ValueError(f"{label}: expected a distinct single-base substitution")
    if not snv and len(ref) == len(alt):
        raise ValueError(f"{label}: non-indel in indel input")
    if snv:
        for column in ("avi_raw", "avi_phred"):
            value = row[column]
            if value is None:
                continue
            try:
                value = float(value)
            except (TypeError, ValueError):
                raise ValueError(f"{label}: invalid {column}") from None
            if math.isnan(value):
                row[column] = None
            elif not math.isfinite(value) or (column == "avi_phred" and value < 0):
                raise ValueError(f"{label}: invalid {column}")
            else:
                row[column] = value
    return row


def audit_variants(rows: list[dict], blocks: dict, *, snv: bool):
    observed, excluded, seen = [], [], {}
    per_block = {key: Counter() for key in blocks}
    columns = SNV_COLUMNS if snv else INDEL_COLUMNS
    for source in rows:
        row = validate_variant(source, blocks, snv=snv)
        key = coordinate_key(row["chrom"], row["block_start"], row["block_end"])
        counts = per_block[key]
        counts["input_rows"] += 1
        frequency = row["AC"] / row["AN"]
        allele = (row["chrom"], row["pos"], row["ref"], row["alt"])
        if allele in seen:
            if any(row[c] != seen[allele][c] for c in columns):
                raise ValueError(f"conflicting duplicate allele: {allele}")
            reason = "duplicate_allele"
        else:
            seen[allele] = row
            reason = "AC_zero" if row["AC"] == 0 else None
        if reason:
            excluded.append({**row, "exclusion_reason": reason})
            counts[reason] += 1
            continue
        row.update(
            block_key=f"{key[0]}:{key[1]}-{key[2]}",
            block_id_avi=blocks[key]["block_id"],
            AF_from_counts=frequency,
            MAF_from_counts=min(frequency, 1 - frequency),
        )
        if snv:
            row["has_avi_score"] = row["avi_phred"] is not None
            counts["observed_scored"] += int(row["has_avi_score"])
            counts["observed_missing_avi"] += int(not row["has_avi_score"])
        counts["observed_alleles"] += 1
        observed.append(row)
    observed.sort(key=lambda r: (r["chrom"], r["pos"], r["ref"], r["alt"]))
    return observed, excluded, per_block


def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snvs", type=Path, default=HERE / "parquet/chr22_small_variants_avi.parquet")
    parser.add_argument("--indels", type=Path, default=HERE / "parquet/chr22_small_variants_indels.parquet")
    parser.add_argument("--block-scores", type=Path, default=HERE / "block_scores.tsv")
    parser.add_argument("--outdir", type=Path, default=HERE / "qc")
    args = parser.parse_args(argv)
    snv_table, indel_table = pq.read_table(args.snvs), pq.read_table(args.indels)
    require_columns(snv_table.column_names, SNV_COLUMNS, str(args.snvs))
    require_columns(indel_table.column_names, INDEL_COLUMNS, str(args.indels))
    if set(snv_table.column_names) & {f.name for f in ADDED_FIELDS}:
        raise ValueError("SNV input already contains QC-derived columns; use the original input")
    blocks = select_blocks(load_block_scores(args.block_scores), pilot_blocks())
    observed, excluded, snv_counts = audit_variants(snv_table.to_pylist(), blocks, snv=True)
    _, _, indel_counts = audit_variants(indel_table.to_pylist(), blocks, snv=False)
    block_rows = []
    for key, block in blocks.items():
        row = {**block, "block_key": f"{key[0]}:{key[1]}-{key[2]}"}
        for prefix, counts in (("snv", snv_counts[key]), ("indel", indel_counts[key])):
            for column in ("input_rows", "duplicate_allele", "AC_zero", "observed_alleles"):
                row[f"{prefix}_{column}"] = counts[column]
        for column in ("observed_scored", "observed_missing_avi"):
            row[f"snv_{column}"] = snv_counts[key][column]
        block_rows.append(row)
    names = ["chr22_observed_snvs.parquet", "chr22_excluded_snvs.tsv", "chr22_block_qc.tsv"]
    inputs = {path.resolve() for path in (args.snvs, args.indels, args.block_scores)}
    if any((args.outdir / name).resolve() in inputs for name in names):
        raise ValueError("output would overwrite an input file")
    schema = snv_table.schema
    for field in ADDED_FIELDS:
        schema = schema.append(field)
    output_table = pa.Table.from_pylist(observed, schema=schema)
    args.outdir.mkdir(parents=True, exist_ok=True)
    pq.write_table(output_table, args.outdir / names[0], compression="zstd")
    write_tsv(args.outdir / names[1], excluded, [*snv_table.column_names, "exclusion_reason"])
    write_tsv(args.outdir / names[2], block_rows, list(block_rows[0]))
    print(f"QC passed: {len(observed):,} observed SNVs; {len(excluded):,} excluded rows; {len(blocks)} blocks. Outputs: {args.outdir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2) from error

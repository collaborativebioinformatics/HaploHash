#!/usr/bin/env python3
"""Rank haploblocks by AVI constraint, for eQTL cross-validation.

Reads block_scores.tsv (static AVI annotation, one row per block) and writes
a ranking: rank 1 = most constrained (highest avi_top1_fraction -- largest
share of the block's possible mutations sit in the AVI top-disruption tier).

avi_top1_fraction is the primary rank key: it's length-normalized, unlike
avi_max (one extreme site can dominate a short block) or avi_top10_mean
(sensitive to n_scored). Ties broken by avi_max.

Usage:
  python3 rank_blocks.py block_scores.tsv block_avi_rank.tsv
"""

from __future__ import annotations

import csv
import sys


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print(f"usage: {argv[0]} block_scores.tsv block_avi_rank.tsv", file=sys.stderr)
        return 1
    in_path, out_path = argv[1], argv[2]

    with open(in_path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    for row in rows:
        row["avi_top1_fraction"] = float(row["avi_top1_fraction"])
        row["avi_max"] = float(row["avi_max"])

    rows.sort(key=lambda r: (-r["avi_top1_fraction"], -r["avi_max"]))

    fieldnames = ["rank", "avi_constraint_percentile", "block_id", "chrom", "start", "end",
                  "length_bp", "avi_top1_fraction", "avi_max", "avi_top10_mean",
                  "avi_top1_count", "avi_top1_per_kb", "avi_top1_mean", "avi_mean", "n_scored"]

    n = len(rows)
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for i, row in enumerate(rows, start=1):
            row["rank"] = i
            row["avi_constraint_percentile"] = round(100 * (n - i + 1) / n, 4)
            writer.writerow({k: row[k] for k in fieldnames})

    print(f"wrote {out_path}: {n} blocks ranked by avi_top1_fraction desc", file=sys.stderr)
    print(f"rank 1 (most constrained): {rows[0]['block_id']} "
          f"top1_fraction={rows[0]['avi_top1_fraction']:.4f}", file=sys.stderr)
    print(f"rank {n} (least constrained): {rows[-1]['block_id']} "
          f"top1_fraction={rows[-1]['avi_top1_fraction']:.4f}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

#!/usr/bin/env python3
"""Test the constraint-vs-eQTL hypothesis: does AVI constraint anti-correlate
with eQTL density?

Joins block_avi_rank.tsv (AVI constraint ranking) with eqtl_by_block.tsv
(GTEx v8 significant eQTLs per block) on block_id, reports Spearman
correlation between avi_top1_fraction and eqtl density, and writes a
rank-vs-eQTL scatter plot.

Usage:
  python3 correlate_avi_eqtl.py block_avi_rank.tsv eqtl_by_block.tsv out_prefix
"""

from __future__ import annotations

import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr


def read_tsv(path: str) -> dict[str, dict]:
    with open(path, newline="") as fh:
        return {row["block_id"]: row for row in csv.DictReader(fh, delimiter="\t")}


def main(argv: list[str]) -> int:
    if len(argv) != 4:
        print(f"usage: {argv[0]} block_avi_rank.tsv eqtl_by_block.tsv out_prefix", file=sys.stderr)
        return 1
    rank_path, eqtl_path, out_prefix = argv[1], argv[2], argv[3]

    rank_rows = read_tsv(rank_path)
    eqtl_rows = read_tsv(eqtl_path)
    joined = [(rank_rows[bid], eqtl_rows[bid]) for bid in rank_rows if bid in eqtl_rows]
    print(f"joined {len(joined)} / {len(rank_rows)} blocks", file=sys.stderr)

    top1_fraction = [float(r["avi_top1_fraction"]) for r, _ in joined]
    avi_mean = [float(r["avi_mean"]) for r, _ in joined]
    avi_rank = [int(r["rank"]) for r, _ in joined]
    eqtl_density = [float(e["eqtl_variants_per_kb"]) for _, e in joined]
    eqtl_count = [int(e["eqtl_variant_count"]) for _, e in joined]

    rho_frac, p_frac = spearmanr(top1_fraction, eqtl_density)
    rho_rank, p_rank = spearmanr(avi_rank, eqtl_count)
    rho_mean, p_mean = spearmanr(avi_mean, eqtl_density)
    print(f"Spearman(avi_top1_fraction, eqtl_variants_per_kb) = {rho_frac:.4f} (p={p_frac:.2e})")
    print(f"Spearman(avi_rank, eqtl_variant_count)            = {rho_rank:.4f} (p={p_rank:.2e})")
    print(f"Spearman(avi_mean [unnormalized], eqtl_variants_per_kb) = {rho_mean:.4f} (p={p_mean:.2e})")

    zero_eqtl = sum(1 for c in eqtl_count if c == 0)
    print(f"blocks with 0 eQTLs: {zero_eqtl} / {len(joined)} ({100*zero_eqtl/len(joined):.1f}%)")

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))

    axes[0].scatter(top1_fraction, eqtl_density, s=4, alpha=0.25, color="#2b6cb0")
    axes[0].set_xlabel("avi_top1_fraction (block constraint, normalized)")
    axes[0].set_ylabel("eQTL variants / kb")
    axes[0].set_title(f"rho={rho_frac:.3f}, p={p_frac:.1e}")

    axes[1].scatter(avi_rank, eqtl_count, s=4, alpha=0.25, color="#c53030")
    axes[1].set_xlabel("AVI constraint rank (1 = most constrained)")
    axes[1].set_ylabel("eQTL variant count")
    axes[1].set_yscale("symlog")
    axes[1].set_title(f"rho={rho_rank:.3f}, p={p_rank:.1e}")

    axes[2].scatter(avi_mean, eqtl_density, s=4, alpha=0.25, color="#2f855a")
    axes[2].set_xlabel("avi_mean (unnormalized, all scored positions)")
    axes[2].set_ylabel("eQTL variants / kb")
    axes[2].set_title(f"rho={rho_mean:.3f}, p={p_mean:.1e}")

    fig.suptitle("AVI constraint vs. GTEx v8 eQTL burden, per haploblock")
    fig.tight_layout()
    fig.savefig(f"{out_prefix}.png", dpi=150)
    print(f"wrote {out_prefix}.png", file=sys.stderr)

    with open(f"{out_prefix}_joined.tsv", "w") as fh:
        fh.write("block_id\tavi_rank\tavi_top1_fraction\tavi_mean\teqtl_variant_count\teqtl_variants_per_kb\n")
        for r, e in joined:
            fh.write(f"{r['block_id']}\t{r['rank']}\t{r['avi_top1_fraction']}\t{r['avi_mean']}\t"
                     f"{e['eqtl_variant_count']}\t{e['eqtl_variants_per_kb']}\n")
    print(f"wrote {out_prefix}_joined.tsv", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

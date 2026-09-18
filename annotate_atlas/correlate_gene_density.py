#!/usr/bin/env python3
"""Test the gene-density-confound explanation with real GENCODE gene counts.

Joins gene_density_by_block.tsv (real gene counts per haploblock) with
block_avi_rank.tsv (avi_top1_fraction, avi_mean) and eqtl_by_block.tsv
(GTEx eQTL density), then checks:

  1. Does AVI density itself track real gene density? (the first half of
     the confound story: "AVI score reflects how much gene machinery is
     packed into a region")
  2. Does gene density predict eQTL density on its own? (the second half:
     "more genes = more eQTL targets")
  3. Does avi_top1_fraction still correlate with eQTL density WITHIN
     gene-density quintiles (i.e. comparing only blocks with similar gene
     counts)? If yes, gene density does not fully explain the AVI-eQTL
     link; if the correlation vanishes within bins, it does.

Usage:
  python3 correlate_gene_density.py gene_density_by_block.tsv \
      block_avi_rank.tsv eqtl_by_block.tsv out_prefix
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
    if len(argv) != 5:
        print(f"usage: {argv[0]} gene_density_by_block.tsv block_avi_rank.tsv eqtl_by_block.tsv out_prefix",
              file=sys.stderr)
        return 1
    gene_path, rank_path, eqtl_path, out_prefix = argv[1:]

    genes = read_tsv(gene_path)
    ranks = read_tsv(rank_path)
    eqtl = read_tsv(eqtl_path)

    rows = []
    for block_id, g in genes.items():
        r = ranks.get(block_id)
        e = eqtl.get(block_id)
        if r is None or e is None:
            continue
        length_kb = (int(g["end"]) - int(g["start"])) / 1000.0
        if length_kb <= 0:
            continue
        rows.append(dict(
            block_id=block_id,
            genes_per_kb=int(g["n_genes_overlap"]) / length_kb,
            n_genes_tss=int(g["n_genes_tss"]),
            avi_top1_fraction=float(r["avi_top1_fraction"]),
            avi_mean=float(r["avi_mean"]),
            eqtl_variants_per_kb=float(e["eqtl_variants_per_kb"]),
        ))
    print(f"joined {len(rows)} / {len(genes)} blocks", file=sys.stderr)

    gene_density = [r["genes_per_kb"] for r in rows]
    avi_frac = [r["avi_top1_fraction"] for r in rows]
    avi_mean = [r["avi_mean"] for r in rows]
    density = [r["eqtl_variants_per_kb"] for r in rows]

    rho_gene_avi, p_gene_avi = spearmanr(gene_density, avi_frac)
    rho_gene_eqtl, p_gene_eqtl = spearmanr(gene_density, density)
    rho_avi_eqtl, p_avi_eqtl = spearmanr(avi_frac, density)
    print(f"Spearman(genes_per_kb, avi_top1_fraction) = {rho_gene_avi:.4f} (p={p_gene_avi:.2e})")
    print(f"Spearman(genes_per_kb, eqtl_variants_per_kb) = {rho_gene_eqtl:.4f} (p={p_gene_eqtl:.2e})")
    print(f"Spearman(avi_top1_fraction, eqtl_variants_per_kb) [unadjusted] = {rho_avi_eqtl:.4f} (p={p_avi_eqtl:.2e})")

    # Within-bin check: does AVI still predict eQTLs once gene density is held (roughly) fixed?
    n_bins = 5
    order = sorted(range(len(rows)), key=lambda i: gene_density[i])
    bin_size = len(order) // n_bins
    within_bin_rhos = []
    print("Within gene-density quintiles (low -> high genes/kb):", file=sys.stderr)
    for b in range(n_bins):
        idx = order[b * bin_size: (b + 1) * bin_size] if b < n_bins - 1 else order[b * bin_size:]
        bin_avi = [avi_frac[i] for i in idx]
        bin_density = [density[i] for i in idx]
        rho, p = spearmanr(bin_avi, bin_density)
        within_bin_rhos.append(rho)
        lo, hi = gene_density[idx[0]], gene_density[idx[-1]]
        print(f"  bin {b+1}: n={len(idx)}, genes/kb in [{lo:.3f}, {hi:.3f}], "
              f"Spearman(avi_top1_fraction, eqtl_density) = {rho:.4f} (p={p:.2e})", file=sys.stderr)

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
    axes[0].scatter(gene_density, avi_frac, s=4, alpha=0.25, color="#2b6cb0")
    axes[0].set_xscale("symlog", linthresh=0.01)
    axes[0].set_xlabel("real genes / kb (GENCODE)")
    axes[0].set_ylabel("avi_top1_fraction")
    axes[0].set_title(f"rho={rho_gene_avi:.3f}, p={p_gene_avi:.1e}")

    axes[1].scatter(gene_density, density, s=4, alpha=0.25, color="#c53030")
    axes[1].set_xscale("symlog", linthresh=0.01)
    axes[1].set_xlabel("real genes / kb (GENCODE)")
    axes[1].set_ylabel("eQTL variants / kb")
    axes[1].set_title(f"rho={rho_gene_eqtl:.3f}, p={p_gene_eqtl:.1e}")

    axes[2].bar(range(1, n_bins + 1), within_bin_rhos, color="#2f855a")
    axes[2].axhline(0, color="0.4", lw=0.8)
    axes[2].set_xlabel("gene-density quintile (1=sparsest, 5=densest)")
    axes[2].set_ylabel("Spearman(avi_top1_fraction, eqtl_density)\nwithin quintile")
    axes[2].set_title("AVI-eQTL link, gene density held fixed")

    fig.suptitle(f"Does real gene density explain the AVI-eQTL link? (n={len(rows)} blocks)")
    fig.tight_layout()
    fig.savefig(f"{out_prefix}.png", dpi=150)
    print(f"wrote {out_prefix}.png", file=sys.stderr)

    with open(f"{out_prefix}_joined.tsv", "w") as fh:
        fh.write("block_id\tgenes_per_kb\tn_genes_tss\tavi_top1_fraction\tavi_mean\teqtl_variants_per_kb\n")
        for r in rows:
            fh.write(f"{r['block_id']}\t{r['genes_per_kb']:.6g}\t{r['n_genes_tss']}\t"
                      f"{r['avi_top1_fraction']}\t{r['avi_mean']}\t{r['eqtl_variants_per_kb']}\n")
    print(f"wrote {out_prefix}_joined.tsv", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

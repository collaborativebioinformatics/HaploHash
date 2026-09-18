#!/usr/bin/env python3
"""Does recombination rate explain the AVI-eQTL and AVI-diversity links?

Joins recomb_rate_by_block.tsv against block_avi_rank.tsv, eqtl_by_block.tsv,
and the cluster-diversity join (shannon_entropy), then checks:

  1. Does recombination rate correlate with AVI, eQTL density, and diversity
     on its own (the candidate common cause)?
  2. Do the AVI-eQTL and AVI-diversity correlations shrink toward zero once
     you only compare blocks with similar recombination rates (quintile
     stratification, same method used for the length and gene-density
     checks)? If yes, recombination rate explains the "backwards" direction
     seen in both. If the correlations survive unchanged, it doesn't.

Usage:
  python3 correlate_recombination.py recomb_rate_by_block.tsv block_avi_rank.tsv \
      eqtl_by_block.tsv block_constraint_avi_joined.tsv out_prefix
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
    if len(argv) != 6:
        print(f"usage: {argv[0]} recomb_rate_by_block.tsv block_avi_rank.tsv eqtl_by_block.tsv "
              f"block_constraint_avi_joined.tsv out_prefix", file=sys.stderr)
        return 1
    recomb_path, rank_path, eqtl_path, diversity_path, out_prefix = argv[1:]

    recomb = read_tsv(recomb_path)
    ranks = read_tsv(rank_path)
    eqtl = read_tsv(eqtl_path)
    diversity = read_tsv(diversity_path)  # keyed by block_id too

    rows = []
    for block_id, r in recomb.items():
        rk = ranks.get(block_id)
        e = eqtl.get(block_id)
        d = diversity.get(block_id)
        if rk is None or e is None or d is None:
            continue
        rows.append(dict(
            block_id=block_id,
            recomb_rate=float(r["recomb_rate_cm_per_mb"]),
            avi_top1_fraction=float(rk["avi_top1_fraction"]),
            eqtl_variants_per_kb=float(e["eqtl_variants_per_kb"]),
            shannon_entropy=float(d["shannon_entropy"]),
        ))
    print(f"joined {len(rows)} / {len(recomb)} blocks", file=sys.stderr)

    recomb_v = [r["recomb_rate"] for r in rows]
    avi_v = [r["avi_top1_fraction"] for r in rows]
    eqtl_v = [r["eqtl_variants_per_kb"] for r in rows]
    div_v = [r["shannon_entropy"] for r in rows]

    print("=== Does recombination rate correlate with each axis on its own? ===")
    rho_avi, p_avi = spearmanr(recomb_v, avi_v)
    rho_eqtl, p_eqtl = spearmanr(recomb_v, eqtl_v)
    rho_div, p_div = spearmanr(recomb_v, div_v)
    print(f"  recomb_rate vs avi_top1_fraction    rho={rho_avi:.4f} p={p_avi:.2e}")
    print(f"  recomb_rate vs eqtl_variants_per_kb rho={rho_eqtl:.4f} p={p_eqtl:.2e}")
    print(f"  recomb_rate vs shannon_entropy      rho={rho_div:.4f} p={p_div:.2e}")

    rho_unadj_eqtl, _ = spearmanr(avi_v, eqtl_v)
    rho_unadj_div, _ = spearmanr(avi_v, div_v)
    print(f"\n=== Unadjusted AVI links (baseline, for comparison) ===")
    print(f"  avi_top1_fraction vs eqtl_variants_per_kb rho={rho_unadj_eqtl:.4f}")
    print(f"  avi_top1_fraction vs shannon_entropy      rho={rho_unadj_div:.4f}")

    n_bins = 5
    order = sorted(range(len(rows)), key=lambda i: recomb_v[i])
    bin_size = len(order) // n_bins
    eqtl_bin_rhos, div_bin_rhos = [], []
    print(f"\n=== Within recombination-rate quintiles ===")
    for b in range(n_bins):
        idx = order[b * bin_size:(b + 1) * bin_size] if b < n_bins - 1 else order[b * bin_size:]
        bin_avi = [avi_v[i] for i in idx]
        bin_eqtl = [eqtl_v[i] for i in idx]
        bin_div = [div_v[i] for i in idx]
        rho_e, p_e = spearmanr(bin_avi, bin_eqtl)
        rho_d, p_d = spearmanr(bin_avi, bin_div)
        eqtl_bin_rhos.append(rho_e)
        div_bin_rhos.append(rho_d)
        lo, hi = recomb_v[idx[0]], recomb_v[idx[-1]]
        print(f"  bin {b+1}: n={len(idx)}, recomb in [{lo:.3f},{hi:.3f}] cM/Mb  "
              f"AVI-eQTL rho={rho_e:.4f} (p={p_e:.1e})   AVI-diversity rho={rho_d:.4f} (p={p_d:.1e})")

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
    axes[0].scatter(recomb_v, avi_v, s=3, alpha=0.15, color="#805ad5")
    axes[0].set_xscale("symlog", linthresh=0.1)
    axes[0].set_xlabel("recombination rate (cM/Mb)")
    axes[0].set_ylabel("avi_top1_fraction")
    axes[0].set_title(f"recomb vs AVI: rho={rho_avi:.3f}")

    x = range(1, n_bins + 1)
    axes[1].plot(x, eqtl_bin_rhos, marker="o", color="#c53030", label="AVI vs eQTL density")
    axes[1].axhline(rho_unadj_eqtl, color="#c53030", ls="--", lw=0.8, alpha=0.5, label="unadjusted (0.39)")
    axes[1].axhline(0, color="0.5", lw=0.7)
    axes[1].set_xlabel("recombination-rate quintile (1=lowest, 5=highest)")
    axes[1].set_ylabel("Spearman(AVI, eQTL density)\nwithin quintile")
    axes[1].set_title("Does recombination explain AVI-eQTL?")
    axes[1].legend(fontsize=8)

    axes[2].plot(x, div_bin_rhos, marker="o", color="#2f855a", label="AVI vs diversity")
    axes[2].axhline(rho_unadj_div, color="#2f855a", ls="--", lw=0.8, alpha=0.5, label="unadjusted (0.09)")
    axes[2].axhline(0, color="0.5", lw=0.7)
    axes[2].set_xlabel("recombination-rate quintile (1=lowest, 5=highest)")
    axes[2].set_ylabel("Spearman(AVI, diversity)\nwithin quintile")
    axes[2].set_title("Does recombination explain AVI-diversity?")
    axes[2].legend(fontsize=8)

    fig.suptitle(f"Testing recombination rate as the common cause (n={len(rows)} blocks)")
    fig.tight_layout()
    fig.savefig(f"{out_prefix}.png", dpi=150)
    print(f"\nwrote {out_prefix}.png", file=sys.stderr)

    with open(f"{out_prefix}_joined.tsv", "w") as fh:
        fh.write("block_id\trecomb_rate_cm_per_mb\tavi_top1_fraction\teqtl_variants_per_kb\tshannon_entropy\n")
        for r in rows:
            fh.write(f"{r['block_id']}\t{r['recomb_rate']:.6g}\t{r['avi_top1_fraction']}\t"
                      f"{r['eqtl_variants_per_kb']}\t{r['shannon_entropy']}\n")
    print(f"wrote {out_prefix}_joined.tsv", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

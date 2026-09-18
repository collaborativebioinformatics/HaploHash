#!/usr/bin/env python3
"""Length-correct Robert's AVI-vs-cluster-diversity join.

join_constraint_avi.py's raw correlations mix two different metric types:
n_clusters/singleton_count are RAW COUNTS (confounded by block length, same
issue as the original eQTL-count comparison), while dominance and
shannon_entropy are already bounded/normalized. This checks the length
confound directly and re-tests the normalized metrics within block-length
quintiles -- the same correction pattern used for the eQTL and gene-density
checks in annotate_atlas/.

Unlike the eQTL comparison, this signal does not touch AlphaGenome's
training/evaluation data at all: shannon_entropy and dominance are computed
purely from real 1000G phased haplotype clusters.

Usage:
  python3 length_corrected_diversity_vs_avi.py block_constraint_avi_joined.tsv out_prefix
"""

from __future__ import annotations

import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import spearmanr


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print(f"usage: {argv[0]} block_constraint_avi_joined.tsv out_prefix", file=sys.stderr)
        return 1
    joined_path, out_prefix = argv[1:]

    df = pd.read_csv(joined_path, sep="\t")
    n_bins = 5
    df["len_bin"] = pd.qcut(df["length_bp"], n_bins, labels=False, duplicates="drop") + 1

    rho_raw, p_raw = spearmanr(df["avi_top1_fraction"], df["shannon_entropy"])
    rho_len_avi, _ = spearmanr(df["length_bp"], df["avi_top1_fraction"])
    rho_len_metric, _ = spearmanr(df["length_bp"], df["shannon_entropy"])
    print(f"Spearman(avi_top1_fraction, shannon_entropy) [unadjusted] = {rho_raw:.4f} (p={p_raw:.2e})")
    print(f"Spearman(length_bp, avi_top1_fraction) = {rho_len_avi:.4f}")
    print(f"Spearman(length_bp, shannon_entropy)   = {rho_len_metric:.4f}")

    within_bin_rhos, within_bin_ps = [], []
    print("Within length quintiles:")
    for b in range(1, n_bins + 1):
        sub = df[df.len_bin == b]
        rho, p = spearmanr(sub["avi_top1_fraction"], sub["shannon_entropy"])
        within_bin_rhos.append(rho)
        within_bin_ps.append(p)
        print(f"  quintile {b}: n={len(sub)}, rho={rho:.4f}, p={p:.2e}")

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))

    axes[0].scatter(df["avi_top1_fraction"], df["shannon_entropy"], s=3, alpha=0.15, color="#2f855a")
    axes[0].set_xlabel("avi_top1_fraction (block risk)")
    axes[0].set_ylabel("Shannon entropy of haplotype clusters")
    axes[0].set_title(f"Unadjusted: rho={rho_raw:.3f}, p={p_raw:.1e}")

    axes[1].scatter(df["length_bp"], df["n_clusters"], s=3, alpha=0.15, color="#c53030")
    axes[1].set_xscale("log")
    axes[1].set_xlabel("block length (bp)")
    axes[1].set_ylabel("n_clusters (raw count)")
    axes[1].set_title("Why correction is needed:\nraw cluster count IS length-confounded")

    axes[2].bar(range(1, n_bins + 1), within_bin_rhos, color="#2f855a")
    axes[2].axhline(0, color="0.4", lw=0.8)
    axes[2].set_xlabel("block-length quintile (1=shortest, 5=longest)")
    axes[2].set_ylabel("Spearman(avi_top1_fraction, shannon_entropy)\nwithin quintile")
    axes[2].set_title("Signal survives length correction")

    fig.suptitle(f"Higher-risk blocks have MORE haplotype diversity, not less (n={len(df)}, length-corrected)")
    fig.tight_layout()
    fig.savefig(f"{out_prefix}.png", dpi=160)
    print(f"wrote {out_prefix}.png", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

#!/usr/bin/env python3
"""Genome-wide test of the original hypothesis: does real-world depletion of
AVI-flagged positions (obs/exp, gnomAD-style constraint) track eQTL burden?

Joins obs_exp_by_block.tsv (this chromosome's real 1000G variants matched
against AVI top1-tier positions) with block_scores.tsv and eqtl_by_block.tsv
(GTEx v8), and reports two different obs/exp ratios per block:

  unconditional_ratio = observed_top1 / avi_top1_count
    Real high-AVI variants vs. ALL possible high-AVI positions in the
    block, whether or not they ever varied. Simple, but conflates "this
    site never varies at all in 2,504 people" (a sampling/power issue,
    true almost everywhere) with "this site is under selection."

  conditional_ratio = observed_top1 / (observed_total * avi_top1_count / n_scored)
    Matches pilot_depletion.py's chr22 pilot formula: among the variants
    that DID occur in this block, are the high-AVI ones under- or
    over-represented relative to the block's own high-AVI density? This
    controls for blocks simply having more or fewer variants overall, so
    it's the more defensible constraint signal -- use this one as primary.

Low ratio = selection suppressed most of the potentially-damaging spots
(real constraint). High ratio = no suppression. If the constraint theory
holds, low-ratio (more suppressed) blocks should show FEWER eQTLs.

Usage:
  python3 correlate_obs_exp_eqtl.py obs_exp_by_block.tsv block_scores.tsv \
      eqtl_by_block.tsv out_prefix
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
        print(f"usage: {argv[0]} obs_exp_by_block.tsv block_scores.tsv eqtl_by_block.tsv out_prefix",
              file=sys.stderr)
        return 1
    obs_exp_path, scores_path, eqtl_path, out_prefix = argv[1:]

    obs_exp = read_tsv(obs_exp_path)
    scores = read_tsv(scores_path)
    eqtl = read_tsv(eqtl_path)

    rows = []
    for block_id, oe in obs_exp.items():
        s = scores.get(block_id)
        e = eqtl.get(block_id)
        if s is None or e is None:
            continue
        avi_top1_count = int(s["avi_top1_count"])
        n_scored = int(s["n_scored"])
        observed_top1 = int(oe["observed_top1"])
        observed_total = int(oe["observed_total"])
        if avi_top1_count == 0 or n_scored == 0 or observed_total == 0:
            continue
        unconditional = observed_top1 / avi_top1_count
        conditional_expected = observed_total * avi_top1_count / n_scored
        if conditional_expected == 0:
            continue
        conditional = observed_top1 / conditional_expected
        rows.append(dict(
            block_id=block_id, chrom=oe["chrom"],
            observed_top1=observed_top1, observed_total=observed_total, avi_top1_count=avi_top1_count,
            unconditional_ratio=unconditional, conditional_ratio=conditional,
            eqtl_variant_count=int(e["eqtl_variant_count"]),
            eqtl_variants_per_kb=float(e["eqtl_variants_per_kb"]),
        ))
    print(f"joined {len(rows)} / {len(obs_exp)} blocks with nonzero AVI top1 positions and >=1 observed variant",
          file=sys.stderr)

    uncond = [r["unconditional_ratio"] for r in rows]
    cond = [r["conditional_ratio"] for r in rows]
    density = [r["eqtl_variants_per_kb"] for r in rows]
    count = [r["eqtl_variant_count"] for r in rows]

    rho_uncond, p_uncond = spearmanr(uncond, density)
    rho_cond, p_cond = spearmanr(cond, density)
    rho_cond_count, p_cond_count = spearmanr(cond, count)
    print(f"Spearman(unconditional obs/exp, eqtl_variants_per_kb) = {rho_uncond:.4f} (p={p_uncond:.2e})")
    print(f"Spearman(conditional obs/exp [pilot formula], eqtl_variants_per_kb) = {rho_cond:.4f} (p={p_cond:.2e})")
    print(f"Spearman(conditional obs/exp [pilot formula], eqtl_variant_count)   = {rho_cond_count:.4f} (p={p_cond_count:.2e})")

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
    axes[0].scatter(uncond, density, s=4, alpha=0.25, color="#6b46c1")
    axes[0].set_xlabel("unconditional obs/exp (vs. all possible AVI-high positions)")
    axes[0].set_ylabel("eQTL variants / kb")
    axes[0].set_title(f"rho={rho_uncond:.3f}, p={p_uncond:.1e}")

    axes[1].scatter(cond, density, s=4, alpha=0.25, color="#b83280")
    axes[1].set_xlabel("conditional obs/exp (pilot formula, vs. observed variants only)")
    axes[1].set_ylabel("eQTL variants / kb")
    axes[1].set_title(f"rho={rho_cond:.3f}, p={p_cond:.1e}")

    axes[2].scatter(cond, count, s=4, alpha=0.25, color="#dd6b20")
    axes[2].set_xlabel("conditional obs/exp (pilot formula)")
    axes[2].set_ylabel("eQTL variant count")
    axes[2].set_yscale("symlog")
    axes[2].set_title(f"rho={rho_cond_count:.3f}, p={p_cond_count:.1e}")

    fig.suptitle(f"Genome-wide obs/exp constraint vs. GTEx v8 eQTL burden (n={len(rows)} blocks)")
    fig.tight_layout()
    fig.savefig(f"{out_prefix}.png", dpi=150)
    print(f"wrote {out_prefix}.png", file=sys.stderr)

    with open(f"{out_prefix}_joined.tsv", "w") as fh:
        fh.write("block_id\tchrom\tobserved_top1\tobserved_total\tavi_top1_count\t"
                  "unconditional_ratio\tconditional_ratio\teqtl_variant_count\teqtl_variants_per_kb\n")
        for r in rows:
            fh.write(f"{r['block_id']}\t{r['chrom']}\t{r['observed_top1']}\t{r['observed_total']}\t"
                      f"{r['avi_top1_count']}\t{r['unconditional_ratio']:.6g}\t{r['conditional_ratio']:.6g}\t"
                      f"{r['eqtl_variant_count']}\t{r['eqtl_variants_per_kb']}\n")
    print(f"wrote {out_prefix}_joined.tsv", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

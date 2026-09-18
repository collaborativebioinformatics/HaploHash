"""Compare observed high-AVI SNVs with a uniform within-block baseline."""

from __future__ import annotations

import argparse
import math
import sys
from collections import Counter
from pathlib import Path

import pyarrow.parquet as pq

from qc_variant_parquet import (
    HERE, audit_variants, coordinate_key, load_block_scores, pilot_blocks,
    select_blocks, write_tsv,
)
from score_blocks import PHRED_THRESHOLD


def summarize(variants: list[dict], blocks: dict) -> tuple[list[dict], dict]:
    observed, excluded, counts = audit_variants(variants, blocks, snv=True)
    if excluded:
        raise ValueError("input contains duplicate or AC=0 rows; run qc_variant_parquet.py first")
    top1 = Counter()
    frequencies = {"high": [], "low": []}
    for row in observed:
        if not row["has_avi_score"]:
            continue
        key = coordinate_key(row["chrom"], row["block_start"], row["block_end"])
        high = row["avi_phred"] >= PHRED_THRESHOLD
        top1[key] += int(high)
        frequencies["high" if high else "low"].append(row["MAF_from_counts"])

    results = []
    for key, block in blocks.items():
        scored = counts[key]["observed_scored"]
        possible, possible_top1 = block["n_scored"], block["avi_top1_count"]
        if top1[key] > possible_top1 or scored - top1[key] > possible - possible_top1:
            raise ValueError(f"observed allele counts exceed possible alleles: {key}")
        expected = scored * possible_top1 / possible if possible else None
        results.append({
            "block_id": block["block_id"], "chrom": key[0], "start": key[1], "end": key[2],
            "n_observed": counts[key]["observed_alleles"],
            "n_observed_scored": scored,
            "n_missing_avi": counts[key]["observed_missing_avi"],
            "n_possible_scored": possible,
            "n_possible_top1": possible_top1,
            "observed_top1": top1[key],
            "expected_top1": expected,
            "observed_expected_ratio": top1[key] / expected if expected else None,
        })
    return results, frequencies


def pooled_totals(rows: list[dict]) -> tuple[int, float, float | None]:
    observed = sum(row["observed_top1"] for row in rows)
    expected = math.fsum(row["expected_top1"] for row in rows if row["expected_top1"] is not None)
    return observed, expected, observed / expected if expected else None


def plot_results(rows: list[dict], frequencies: dict, outdir: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, PercentFormatter

    plt.style.use(HERE / "paper.mplstyle")
    blue, orange = plt.rcParams["axes.prop_cycle"].by_key()["color"][:2]
    gray = "#B0B0B0"
    observed, expected, ratio = pooled_totals(rows)
    ratio_text = f"{ratio:.2f}" if ratio is not None else "undefined"
    fig, ax = plt.subplots(figsize=(6, 3.6), layout="constrained")
    for i, row in enumerate(rows):
        if row["expected_top1"] is not None:
            ax.plot([row["observed_top1"], row["expected_top1"]], [i, i], color=gray, lw=0.8)
    ax.scatter([r["observed_top1"] for r in rows], range(len(rows)), color=blue, s=14,
               label=f"Observed (total {observed:,})", zorder=3)
    defined = [(i, r["expected_top1"]) for i, r in enumerate(rows) if r["expected_top1"] is not None]
    ax.scatter([v for _, v in defined], [i for i, _ in defined], color=orange, marker="|", s=55,
               linewidths=1.1, label=f"Expected (total {expected:.1f})", zorder=3)
    ax.set_yticks(range(len(rows)), [f"{r['start']:,}–{r['end']:,}" for r in rows])
    ax.invert_yaxis()
    ax.set_ylabel("Chr22 block (0-based interval)")
    ax.set_xlabel(f"Alternate alleles with AVI PHRED ≥ {PHRED_THRESHOLD:g}")
    ax.set_xlim(left=-0.8)
    ax.grid(axis="y", visible=False)
    ax.tick_params(axis="y", length=0)
    ax.legend(loc="lower right", title=f"Pooled O/E = {ratio_text}")
    ax.set_title("Observed and expected high-AVI variants")
    fig.savefig(outdir / "observed_vs_expected.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(3.5, 2.7), layout="constrained")
    all_frequencies = frequencies["high"] + frequencies["low"]
    positive = [v for v in all_frequencies if v > 0]
    minimum = min(positive, default=1e-4)
    has_zero = any(v == 0 for v in all_frequencies)
    left = 0 if has_zero else minimum * 0.7
    for category, label, color in (("high", "AVI ≥ 20", orange), ("low", "AVI < 20", blue)):
        values = sorted(frequencies[category])
        if values:
            cumulative = [i / len(values) for i in range(1, len(values) + 1)]
            ax.step([left, *values, 0.5], [0, *cumulative, 1], where="post",
                    color=color, lw=1.1, label=f"{label} (n = {len(values):,})")
    if has_zero:
        ax.set_xscale("symlog", linthresh=minimum / 2)
    else:
        ax.set_xscale("log")
    ax.set_xlim(left, 0.5)
    ax.set_ylim(0, 1.02)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{100 * value:g}%"))
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
    ax.set_xlabel("Minor allele frequency")
    ax.set_ylabel("Cumulative fraction of SNVs")
    ax.set_title("Chr22 allele-frequency distributions")
    if all_frequencies:
        ax.legend(loc="lower right", frameon=False)
    else:
        ax.text(0.5, 0.5, "No scored observed SNVs", transform=ax.transAxes, ha="center")
    fig.savefig(outdir / "allele_frequency_by_avi.png")
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snvs", type=Path, default=HERE / "qc/chr22_observed_snvs.parquet")
    parser.add_argument("--block-scores", type=Path, default=HERE / "qc/chr22_block_qc.tsv")
    parser.add_argument("--outdir", type=Path, default=HERE / "depletion")
    args = parser.parse_args(argv)
    blocks = select_blocks(load_block_scores(args.block_scores), pilot_blocks())
    rows, frequencies = summarize(pq.read_table(args.snvs).to_pylist(), blocks)
    names = ["chr22_depletion.tsv", "observed_vs_expected.png", "allele_frequency_by_avi.png"]
    if any((args.outdir / name).resolve() in {args.snvs.resolve(), args.block_scores.resolve()} for name in names):
        raise ValueError("output would overwrite an input file")
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.outdir / names[0], rows, list(rows[0]))
    plot_results(rows, frequencies, args.outdir)
    observed, expected, ratio = pooled_totals(rows)
    ratio_text = f"{ratio:.3f}" if ratio is not None else "undefined"
    missing = sum(row["n_missing_avi"] for row in rows)
    print(f"{observed:,} observed / {expected:.2f} expected = {ratio_text}; {missing} unscored SNVs omitted. Outputs: {args.outdir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2) from error

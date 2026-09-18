#!/usr/bin/env python3
"""Join haploblock constraint statistics to block-level AlphaGenome Atlas AVI scores,
write a join report, and plot a 3 x 6 matrix of descriptive scatterplots.

Inputs (TSV):
  block_avi_scores.tsv     block_id, chrom, start, end, length_bp, ..., avi_* columns
  block_cluster_stats.tsv  chr, block, start, end, block_length, n_clusters, ...

The two files carry incompatible block-name strings:
    block_avi_scores.block_id      chr1_1583825_1958353     (underscore before end)
    block_cluster_stats.block      chr1_1583825-1958353     (hyphen before end)
so the join is performed on the canonical coordinate key (chrom, start, end),
which both files carry as separate columns. A raw-string join is attempted first
and its failure is recorded in the report.

Usage:
  python3 join_constraint_avi.py [--indir DIR] [--outdir DIR]

Outputs, written to --outdir:
  block_constraint_avi_joined.tsv       merged table
  join_report.txt                       join diagnostics and per-panel statistics
  scatter_matrix_3x6.png / .pdf         large exploratory matrix
  fig_constraint_vs_avi_double.pdf/.png Nature double-column (180 mm) figure
"""

from __future__ import annotations

import argparse
import datetime as _dt
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, MaxNLocator, SymmetricalLogLocator
from scipy import stats

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #

AVI_FILE = "block_avi_scores.tsv"
CLUSTER_FILE = "block_cluster_stats.tsv"

# X axis: block constraint measures (rows of the matrix)
X_MEASURES = [
    ("n_clusters", "Cluster count", "log"),
    ("singleton_count", "Singleton clusters", "symlog"),
    ("dominance", "Dominance", "linear"),
]

# Y axis: block AVI measures (columns of the matrix)
Y_MEASURES = [
    ("avi_max", "avi_max", "linear"),
    ("avi_top10_mean", "avi_top10_mean", "linear"),
    ("avi_top1_count", "avi_top1_count", "symlog"),
    ("avi_top1_fraction", "avi_top1_fraction", "symlog"),
    ("avi_top1_per_kb", "avi_top1_per_kb", "symlog"),
    ("avi_top1_mean", "avi_top1_mean", "linear"),
]

SIZE_COL = "block_length"

# symlog linear-threshold per column (region around zero rendered linearly)
LINTHRESH = {
    "singleton_count": 1.0,
    "avi_top1_count": 1.0,
    "avi_top1_fraction": 1e-4,
    "avi_top1_per_kb": 1e-1,
}

MARK = "#2E5C8A"      # single-series mark colour
TREND = "#C1462C"     # binned-median trend
INK = "#1A1A1A"
MUTED = "#6B6B6B"
GRID = "#DDDDDD"

PANEL_LETTERS = "abcdefghijklmnopqrstuvwxyz"


# --------------------------------------------------------------------------- #
# Join
# --------------------------------------------------------------------------- #

def canonical_key(chrom: pd.Series, start: pd.Series, end: pd.Series) -> pd.Series:
    return (
        chrom.astype(str).str.strip()
        + ":"
        + start.astype("int64").astype(str)
        + "-"
        + end.astype("int64").astype(str)
    )


def load_inputs(indir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    avi = pd.read_csv(indir / AVI_FILE, sep="\t")
    clust = pd.read_csv(indir / CLUSTER_FILE, sep="\t")
    avi["_key"] = canonical_key(avi["chrom"], avi["start"], avi["end"])
    clust["_key"] = canonical_key(clust["chr"], clust["start"], clust["end"])
    return avi, clust


def do_join(avi: pd.DataFrame, clust: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Outer-join on the canonical key; return merged frame and diagnostics."""
    diag: dict = {}

    diag["n_rows_avi"] = len(avi)
    diag["n_rows_cluster"] = len(clust)

    # Raw-string join attempt, recorded for the report.
    raw_avi = set(avi["block_id"].astype(str))
    raw_clust = set(clust["block"].astype(str))
    diag["raw_string_overlap"] = len(raw_avi & raw_clust)
    diag["raw_example_avi"] = sorted(raw_avi)[0] if raw_avi else ""
    diag["raw_example_cluster"] = sorted(raw_clust)[0] if raw_clust else ""

    # Duplicate keys
    diag["dup_keys_avi"] = int(avi["_key"].duplicated().sum())
    diag["dup_keys_cluster"] = int(clust["_key"].duplicated().sum())

    merged = avi.merge(
        clust.drop(columns=["start", "end"]),
        on="_key",
        how="outer",
        indicator=True,
        suffixes=("_avi", "_cl"),
    )

    diag["n_both"] = int((merged["_merge"] == "both").sum())
    diag["n_avi_only"] = int((merged["_merge"] == "left_only").sum())
    diag["n_cluster_only"] = int((merged["_merge"] == "right_only").sum())

    diag["avi_only_keys"] = merged.loc[merged["_merge"] == "left_only", "_key"].tolist()
    diag["cluster_only_keys"] = merged.loc[
        merged["_merge"] == "right_only", "_key"
    ].tolist()

    inner = merged[merged["_merge"] == "both"].copy()

    # Coordinate-length cross-check between the two sources.
    both_len = inner["length_bp"].notna() & inner["block_length"].notna()
    diag["length_checked"] = int(both_len.sum())
    diag["length_agree"] = int(
        (inner.loc[both_len, "length_bp"] == inner.loc[both_len, "block_length"]).sum()
    )
    diag["length_disagree_examples"] = inner.loc[
        both_len
        & (inner["length_bp"] != inner["block_length"]),
        ["_key", "length_bp", "block_length"],
    ].head(10)

    # Per-chromosome match counts
    chrom_all = merged["chrom"].fillna(merged["chr"])
    diag["per_chrom"] = (
        pd.DataFrame({"chrom": chrom_all, "status": merged["_merge"].astype(str)})
        .pivot_table(index="chrom", columns="status", aggfunc=len, fill_value=0)
        .rename(columns={"both": "matched", "left_only": "avi_only",
                         "right_only": "cluster_only"})
    )

    inner = inner.drop(columns=["_merge"])
    return inner, diag


# --------------------------------------------------------------------------- #
# Plot helpers
# --------------------------------------------------------------------------- #

def marker_sizes(length: np.ndarray, smin: float, smax: float) -> np.ndarray:
    """Map block length onto marker area, log-spaced so four orders of magnitude fit."""
    lg = np.log10(np.clip(length.astype(float), 1.0, None))
    lo, hi = np.nanpercentile(lg, 0.5), np.nanpercentile(lg, 99.5)
    if hi <= lo:
        return np.full(lg.shape, smin)
    t = np.clip((lg - lo) / (hi - lo), 0.0, 1.0)
    return smin + t * (smax - smin)


def apply_scale(ax, axis: str, scale: str, col: str) -> None:
    if scale == "log":
        setter = ax.set_xscale if axis == "x" else ax.set_yscale
        setter("log")
    elif scale == "symlog":
        lt = LINTHRESH.get(col, 1.0)
        if axis == "x":
            ax.set_xscale("symlog", linthresh=lt, linscale=0.35)
        else:
            ax.set_yscale("symlog", linthresh=lt, linscale=0.35)


def thin_ticks(ax, axis: str, scale: str, col: str, nbins: int,
               minor: bool = True, log_numticks: int = 5) -> None:
    target = ax.xaxis if axis == "x" else ax.yaxis
    if scale == "log":
        target.set_major_locator(LogLocator(base=10, numticks=log_numticks))
        if minor:
            target.set_minor_locator(
                LogLocator(base=10, subs="auto", numticks=log_numticks)
            )
        else:
            target.set_minor_locator(matplotlib.ticker.NullLocator())
        target.set_minor_formatter(matplotlib.ticker.NullFormatter())
    elif scale == "symlog":
        # SymmetricalLogLocator puts a tick at 0 and another at linthresh; at
        # small panel sizes those two labels collide (rendering as "010^0").
        # Place 0 explicitly, then start the decades one power above linthresh.
        lt = LINTHRESH.get(col, 1.0)
        lo, hi = (ax.get_xlim() if axis == "x" else ax.get_ylim())
        k0 = int(np.ceil(np.log10(lt))) + 1
        k1 = int(np.floor(np.log10(max(hi, lt * 10))))
        ticks = [0.0] + [10.0 ** k for k in range(k0, k1 + 1)]
        target.set_major_locator(matplotlib.ticker.FixedLocator(ticks))
        target.set_minor_locator(matplotlib.ticker.NullLocator())
        target.set_minor_formatter(matplotlib.ticker.NullFormatter())
    else:
        target.set_major_locator(MaxNLocator(nbins=nbins, prune="both"))


def binned_median(x: np.ndarray, y: np.ndarray, scale: str, nbins: int = 20):
    """Median of y within equal-count bins of x; returns bin centres and medians."""
    if len(x) < nbins * 5:
        return None, None
    qs = np.linspace(0.0, 1.0, nbins + 1)
    edges = np.unique(np.nanquantile(x, qs))
    if len(edges) < 3:
        return None, None
    idx = np.clip(np.digitize(x, edges[1:-1]), 0, len(edges) - 2)
    cx, cy = [], []
    for b in range(len(edges) - 1):
        sel = idx == b
        if sel.sum() < 10:
            continue
        cx.append(np.nanmedian(x[sel]))
        cy.append(np.nanmedian(y[sel]))
    if len(cx) < 3:
        return None, None
    return np.asarray(cx), np.asarray(cy)


def panel(ax, df, xcol, ycol, xscale, yscale, smin, smax, dot_alpha,
          rho_fs, trend_lw, show_trend=True):
    """Draw one scatter panel; return (n, spearman rho, p)."""
    sub = df[[xcol, ycol, SIZE_COL]].replace([np.inf, -np.inf], np.nan).dropna()
    if xscale == "log":
        sub = sub[sub[xcol] > 0]
    if yscale == "log":
        sub = sub[sub[ycol] > 0]

    x = sub[xcol].to_numpy(float)
    y = sub[ycol].to_numpy(float)
    s = marker_sizes(sub[SIZE_COL].to_numpy(), smin, smax)

    order = np.argsort(s)  # small blocks first, large drawn on top and visible
    ax.scatter(
        x[order], y[order], s=s[order], c=MARK, alpha=dot_alpha,
        linewidths=0, rasterized=True, zorder=2,
    )

    if show_trend:
        bx, by = binned_median(x, y, xscale)
        if bx is not None:
            ax.plot(bx, by, color=TREND, lw=trend_lw, zorder=4,
                    solid_capstyle="round")
            ax.plot(bx, by, color="white", lw=trend_lw + 1.0, zorder=3,
                    solid_capstyle="round")

    apply_scale(ax, "x", xscale, xcol)
    apply_scale(ax, "y", yscale, ycol)

    # Symlog renders a mirrored negative branch by default. Both constraint and
    # AVI measures are non-negative, so clamp the lower bound to zero; otherwise
    # the data are squeezed into the top of an axis that is mostly empty.
    if xscale == "symlog" and np.nanmin(x) >= 0:
        ax.set_xlim(left=0, right=np.nanmax(x) * 1.6)
    if yscale == "symlog" and np.nanmin(y) >= 0:
        ax.set_ylim(bottom=0, top=np.nanmax(y) * 1.6)
    # Log axes otherwise clip the lowest stripe of points against the spine.
    if xscale == "log":
        ax.set_xlim(left=np.nanmin(x) * 0.65, right=np.nanmax(x) * 1.55)
    if yscale == "log":
        ax.set_ylim(bottom=np.nanmin(y) * 0.65, top=np.nanmax(y) * 1.55)

    rho, p = (np.nan, np.nan)
    if len(x) > 2:
        rho, p = stats.spearmanr(x, y)

    ax.grid(True, color=GRID, lw=0.4, alpha=0.9, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
        ax.spines[side].set_linewidth(0.6)

    if np.isfinite(rho):
        ax.text(
            0.04, 0.955, f"$\\rho$ = {rho:+.2f}", transform=ax.transAxes,
            ha="left", va="top", fontsize=rho_fs, color=INK,
            bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.75),
            zorder=5,
        )
    return len(x), rho, p


def size_legend_handles(df, smin, smax, fs):
    lengths = [1e4, 1e5, 1e6]
    lg = np.log10(np.clip(df[SIZE_COL].to_numpy(float), 1, None))
    lo, hi = np.nanpercentile(lg, 0.5), np.nanpercentile(lg, 99.5)
    handles, labels = [], []
    for L in lengths:
        t = np.clip((np.log10(L) - lo) / (hi - lo), 0, 1)
        area = smin + t * (smax - smin)
        handles.append(
            Line2D([], [], marker="o", linestyle="none", color=MARK,
                   markersize=np.sqrt(area), alpha=0.65, markeredgewidth=0)
        )
        labels.append(f"{L/1e3:,.0f} kb" if L < 1e6 else f"{L/1e6:,.0f} Mb")
    handles.append(Line2D([], [], color=TREND, lw=1.4))
    labels.append("binned median")
    return handles, labels


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #

def build_matrix(df, outpath_stem: Path, *, width_in, height_in, base_fs,
                 tick_fs, header_fs, rho_fs, smin, smax, dot_alpha, trend_lw,
                 nbins_ticks, letters, dpi, wspace, hspace, minor_ticks,
                 row_label_pad, legend_y, log_numticks, margins):
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": base_fs,
        "axes.linewidth": 0.6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.dpi": dpi,
    })

    nrow, ncol = len(X_MEASURES), len(Y_MEASURES)
    # X measure is constant along a row, Y measure constant down a column, so
    # sharing on those axes makes every column directly comparable.
    fig, axes = plt.subplots(
        nrow, ncol, figsize=(width_in, height_in), squeeze=False,
        sharex="row", sharey="col",
        gridspec_kw={"wspace": wspace, "hspace": hspace, **margins},
    )

    stats_rows = []
    k = 0
    for i, (xcol, xlabel, xscale) in enumerate(X_MEASURES):
        for j, (ycol, ylabel, yscale) in enumerate(Y_MEASURES):
            ax = axes[i][j]
            n, rho, p = panel(
                ax, df, xcol, ycol, xscale, yscale,
                smin, smax, dot_alpha, rho_fs, trend_lw,
            )
            thin_ticks(ax, "x", xscale, xcol, nbins_ticks, minor_ticks,
                       log_numticks)
            thin_ticks(ax, "y", yscale, ycol, nbins_ticks, minor_ticks,
                       log_numticks)
            ax.tick_params(axis="both", labelsize=tick_fs, length=2.2,
                           width=0.5, pad=1.6, colors=MUTED,
                           labelbottom=True, labelleft=True)
            for lbl in ax.get_xticklabels() + ax.get_yticklabels():
                lbl.set_color(INK)

            if letters:
                ax.set_title(
                    PANEL_LETTERS[k], loc="left", fontsize=header_fs,
                    fontweight="bold", color=INK, pad=2.0,
                )
            if i == 0:
                ax.annotate(
                    ylabel, xy=(0.5, 1.0), xytext=(0, 14 if letters else 6),
                    xycoords="axes fraction", textcoords="offset points",
                    ha="center", va="bottom", fontsize=header_fs,
                    fontweight="bold", color=INK,
                )
            if j == 0:
                ax.annotate(
                    xlabel, xy=(0, 0.5), xytext=(-row_label_pad, 0),
                    xycoords="axes fraction", textcoords="offset points",
                    ha="center", va="center", rotation=90,
                    fontsize=header_fs, fontweight="bold", color=INK,
                )
            stats_rows.append({
                "panel": PANEL_LETTERS[k],
                "x": xcol, "y": ycol, "n": n,
                "spearman_rho": rho, "p_value": p,
            })
            k += 1

    handles, labels = size_legend_handles(df, smin, smax, tick_fs)
    leg = fig.legend(
        handles, labels, loc="lower center", ncol=4, frameon=False,
        fontsize=tick_fs + 1.0, handletextpad=0.6, columnspacing=2.0,
        borderaxespad=0.0, bbox_to_anchor=(0.5, legend_y),
        title="block length encoded as dot area", title_fontsize=tick_fs + 1.0,
    )
    leg.get_title().set_color(INK)

    for ext in ("png", "pdf"):
        fig.savefig(f"{outpath_stem}.{ext}", dpi=dpi, facecolor="white")
    plt.close(fig)
    return pd.DataFrame(stats_rows)


# --------------------------------------------------------------------------- #
# Report
# --------------------------------------------------------------------------- #

def write_report(path: Path, indir: Path, diag: dict, merged: pd.DataFrame,
                 stats_df: pd.DataFrame) -> None:
    L = []
    w = L.append
    w("=" * 78)
    w("JOIN REPORT: block constraint statistics x block AVI scores")
    w("=" * 78)
    w(f"generated          : {_dt.datetime.now().isoformat(timespec='seconds')}")
    w(f"input directory    : {indir}")
    w(f"left  (AVI)        : {AVI_FILE}   rows = {diag['n_rows_avi']:,}")
    w(f"right (constraint) : {CLUSTER_FILE}   rows = {diag['n_rows_cluster']:,}")
    w("")

    w("-" * 78)
    w("1. KEY COMPATIBILITY")
    w("-" * 78)
    w("A raw join on the block-name strings was attempted first and FAILED:")
    w(f"  block_avi_scores.block_id     example: {diag['raw_example_avi']}")
    w(f"  block_cluster_stats.block     example: {diag['raw_example_cluster']}")
    w(f"  shared raw strings            : {diag['raw_string_overlap']:,}")
    w("  cause: the AVI file separates start and end with '_', the constraint")
    w("         file with '-'. The strings are therefore never equal.")
    w("")
    w("The join was performed on the canonical coordinate key chrom:start-end,")
    w("which both files carry as separate columns.")
    w(f"  duplicate keys, AVI file          : {diag['dup_keys_avi']:,}")
    w(f"  duplicate keys, constraint file   : {diag['dup_keys_cluster']:,}")
    w("")

    w("-" * 78)
    w("2. JOIN RESULT")
    w("-" * 78)
    tot_l, tot_r = diag["n_rows_avi"], diag["n_rows_cluster"]
    w(f"  matched (both files)              : {diag['n_both']:,}"
      f"   ({100*diag['n_both']/tot_l:.2f}% of AVI,"
      f" {100*diag['n_both']/tot_r:.2f}% of constraint)")
    w(f"  AVI only, no constraint record    : {diag['n_avi_only']:,}")
    w(f"  constraint only, no AVI record    : {diag['n_cluster_only']:,}")
    w("")
    if diag["avi_only_keys"]:
        w("  unmatched AVI blocks:")
        for k in diag["avi_only_keys"][:50]:
            w(f"    {k}")
        if len(diag["avi_only_keys"]) > 50:
            w(f"    ... and {len(diag['avi_only_keys'])-50:,} more")
    else:
        w("  unmatched AVI blocks: none")
    w("")
    if diag["cluster_only_keys"]:
        w("  unmatched constraint blocks:")
        for k in diag["cluster_only_keys"][:50]:
            w(f"    {k}")
        if len(diag["cluster_only_keys"]) > 50:
            w(f"    ... and {len(diag['cluster_only_keys'])-50:,} more")
    else:
        w("  unmatched constraint blocks: none")
    w("")

    w("-" * 78)
    w("3. CONSISTENCY CHECK ON MATCHED ROWS")
    w("-" * 78)
    w(f"  rows with block length in both    : {diag['length_checked']:,}")
    w(f"  length_bp == block_length         : {diag['length_agree']:,}")
    w(f"  disagreements                     : "
      f"{diag['length_checked'] - diag['length_agree']:,}")
    if len(diag["length_disagree_examples"]):
        w("  first disagreements:")
        for _, r in diag["length_disagree_examples"].iterrows():
            w(f"    {r['_key']}  avi={r['length_bp']}  constraint={r['block_length']}")
    w("")

    w("-" * 78)
    w("4. MATCHES PER CHROMOSOME")
    w("-" * 78)
    pc = diag["per_chrom"]
    for col in ("matched", "avi_only", "cluster_only"):
        if col not in pc.columns:
            pc[col] = 0
    pc = pc[["matched", "avi_only", "cluster_only"]]
    order = sorted(
        pc.index,
        key=lambda c: (99 if not c[3:].isdigit() else int(c[3:]), c),
    )
    w(f"  {'chrom':<8}{'matched':>10}{'avi_only':>10}{'cluster_only':>14}")
    for c in order:
        r = pc.loc[c]
        w(f"  {c:<8}{int(r['matched']):>10,}{int(r['avi_only']):>10,}"
          f"{int(r['cluster_only']):>14,}")
    w(f"  {'TOTAL':<8}{int(pc['matched'].sum()):>10,}"
      f"{int(pc['avi_only'].sum()):>10,}{int(pc['cluster_only'].sum()):>14,}")
    w("")

    w("-" * 78)
    w("5. MISSING VALUES IN PLOTTED MEASURES (matched rows)")
    w("-" * 78)
    w(f"  {'column':<22}{'n_missing':>12}{'n_zero':>10}{'min':>16}{'max':>16}")
    for col, _, _ in X_MEASURES + Y_MEASURES + [(SIZE_COL, "", "")]:
        s = merged[col]
        w(f"  {col:<22}{int(s.isna().sum()):>12,}{int((s == 0).sum()):>10,}"
          f"{s.min():>16,.6g}{s.max():>16,.6g}")
    w("")
    w("  note: rows with a missing value in either axis are dropped from that")
    w("        panel only; the per-panel n is given in section 6.")
    w("")

    w("-" * 78)
    w("6. PER-PANEL SAMPLE SIZE AND RANK CORRELATION")
    w("-" * 78)
    w(f"  {'panel':<7}{'x':<18}{'y':<20}{'n':>9}{'spearman':>11}{'p':>12}")
    for _, r in stats_df.iterrows():
        p = r["p_value"]
        pstr = "< 1e-300" if (pd.notna(p) and p == 0) else f"{p:.3g}"
        w(f"  {r['panel']:<7}{r['x']:<18}{r['y']:<20}{int(r['n']):>9,}"
          f"{r['spearman_rho']:>11.3f}{pstr:>12}")
    w("")
    w("=" * 78)
    w("END OF REPORT")
    w("=" * 78)

    path.write_text("\n".join(L) + "\n", encoding="utf-8")


# --------------------------------------------------------------------------- #

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", default=".", type=Path)
    ap.add_argument("--outdir", default=".", type=Path)
    args = ap.parse_args(argv)

    indir, outdir = args.indir.resolve(), args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    avi, clust = load_inputs(indir)
    merged, diag = do_join(avi, clust)

    merged_out = merged.drop(columns=["_key"]).copy()
    merged_out.insert(0, "block_key", merged["_key"].values)
    merged_out.to_csv(outdir / "block_constraint_avi_joined.tsv",
                      sep="\t", index=False)

    # Large exploratory matrix (screen reading).
    stats_df = build_matrix(
        merged, outdir / "scatter_matrix_3x6",
        width_in=17.0, height_in=9.2, base_fs=9, tick_fs=7.5, header_fs=10,
        rho_fs=8, smin=0.8, smax=20.0, dot_alpha=0.07, trend_lw=1.7,
        nbins_ticks=5, letters=False, dpi=200,
        wspace=0.30, hspace=0.30, minor_ticks=True, row_label_pad=46,
        legend_y=0.012, log_numticks=5,
        margins={"left": 0.055, "right": 0.995, "top": 0.935, "bottom": 0.115},
    )

    # Nature double-column figure: 180 mm wide.
    build_matrix(
        merged, outdir / "fig_constraint_vs_avi_double",
        width_in=180 / 25.4, height_in=118 / 25.4, base_fs=6, tick_fs=4.8,
        header_fs=6, rho_fs=4.8, smin=0.35, smax=7.0, dot_alpha=0.06,
        trend_lw=0.9, nbins_ticks=5, letters=True, dpi=600,
        wspace=0.52, hspace=0.42, minor_ticks=False, row_label_pad=26,
        legend_y=0.010, log_numticks=4,
        margins={"left": 0.072, "right": 0.992, "top": 0.915, "bottom": 0.125},
    )

    write_report(outdir / "join_report.txt", indir, diag, merged, stats_df)

    print(f"matched {diag['n_both']:,} blocks "
          f"({diag['n_avi_only']:,} AVI-only, {diag['n_cluster_only']:,} constraint-only)")
    print(f"wrote outputs to {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Cluster count versus AVI block scores: ordinary least squares.

x = n_clusters          (block constraint, from block_cluster_stats.tsv)
y = avi_top1_per_kb     (block AVI density, from block_avi_scores.tsv)
y = avi_mean            (block mean AVI, from block_avi_scores.tsv)

Every response in Y_TARGETS is modelled separately against the same predictor,
with its own figures; adding another AVI column means adding it to that list.

Three model variants are fitted per response:
  A  raw y on raw x                 y ~ n_clusters
  B  raw y on log10 x               y ~ log10(n_clusters)
  C  log10 y on log10 x             log10(y) ~ log10(n_clusters)

Variant C is the log-log model and its slope is an elasticity. Blocks with
y == 0 have no logarithm and are dropped from variant C only.

Inputs are the same two TSVs used by join_constraint_avi.py; the join is redone
here on the canonical coordinate key so this script stands alone.

Usage:
  python3 regress_clusters_vs_avi_per_kb.py [--indir DIR] [--outdir DIR]

Outputs, written to --outdir:
  regression_report.txt                   coefficients, CIs, R2
  fitted_values.tsv                       fitted lines with 95% confidence bands
  fig_linear_<y>.png                      OLS fit, 180 mm wide
  panel_linear_<y>_a_raw_x.png            variant A, full y range
  panel_linear_<y>_b_log10_x.png          variant B, full y range
  panel_linear_<y>_c_loglog.png           variant C, log-log
  panel_linear_<y>_d_raw_x_zoom.png       variant A, y truncated at the 99th pct
  panel_linear_<y>_e_log10_x_zoom.png     variant B, y truncated at the 99th pct
  where <y> is the response column. The zoomed panels appear only for a response
  whose maximum exceeds its 99th percentile by more than ZOOM_RATIO.
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
import statsmodels.api as sm
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, MaxNLocator
from scipy import stats

AVI_FILE = "block_avi_scores.tsv"
CLUSTER_FILE = "block_cluster_stats.tsv"

XCOL, SIZE_COL = "n_clusters", "block_length"
# Every AVI response modelled against the same predictor.
Y_TARGETS = ["avi_top1_per_kb", "avi_mean"]
# A zoomed companion panel is drawn only where the top of the y range
# dwarfs the bulk of the data; below this ratio it adds nothing.
ZOOM_RATIO = 2.0
# Output image formats. Add "pdf" here for vector copies.
FORMATS = ("png",)

MARK = "#2E5C8A"      # points
LINFIT = "#C1462C"    # fitted line
INK = "#1A1A1A"
MUTED = "#6B6B6B"
GRID = "#DDDDDD"


# --------------------------------------------------------------------------- #
# Data
# --------------------------------------------------------------------------- #

def load_and_join(indir: Path) -> pd.DataFrame:
    avi = pd.read_csv(indir / AVI_FILE, sep="\t")
    clust = pd.read_csv(indir / CLUSTER_FILE, sep="\t")
    key = lambda c, s, e: (c.astype(str).str.strip() + ":"
                           + s.astype("int64").astype(str) + "-"
                           + e.astype("int64").astype(str))
    avi["_key"] = key(avi["chrom"], avi["start"], avi["end"])
    clust["_key"] = key(clust["chr"], clust["start"], clust["end"])
    return avi.merge(clust.drop(columns=["start", "end"]), on="_key", how="inner")


# --------------------------------------------------------------------------- #
# Fitting
# --------------------------------------------------------------------------- #

def fit_variant(x: np.ndarray, y: np.ndarray, xgrid: np.ndarray,
                log_y: bool = False) -> dict:
    """OLS of y on x, with classical and HC3 covariance and a fitted line."""
    X = sm.add_constant(x[:, None])
    lin = sm.OLS(y, X).fit()
    # The residual variance rises steeply with x, so the classical standard
    # errors are optimistic; HC3 is the column to quote.
    lin_hc3 = lin.get_robustcov_results(cov_type="HC3")

    pred = lin.get_prediction(sm.add_constant(xgrid[:, None])).summary_frame(
        alpha=0.05)
    fit = pred["mean"].to_numpy()
    lo = pred["mean_ci_lower"].to_numpy()
    hi = pred["mean_ci_upper"].to_numpy()
    # When the response is logged the line is fitted in log10 space; the plot
    # shows it back-transformed on a log y axis, where it is still straight.
    back = (lambda v: 10.0 ** v) if log_y else (lambda v: v)
    return {
        "lin": lin, "lin_hc3": lin_hc3,
        "xgrid": xgrid,
        "lin_fit": fit, "lin_lo": lo, "lin_hi": hi,
        "fit_plot": back(fit), "lo_plot": back(lo), "hi_plot": back(hi),
        "log_y": log_y,
        "n": len(x),
    }


# --------------------------------------------------------------------------- #
# Plotting
# --------------------------------------------------------------------------- #

def marker_sizes(length: np.ndarray, smin: float, smax: float) -> np.ndarray:
    lg = np.log10(np.clip(length.astype(float), 1.0, None))
    lo, hi = np.nanpercentile(lg, 0.5), np.nanpercentile(lg, 99.5)
    t = np.clip((lg - lo) / (hi - lo), 0.0, 1.0)
    return smin + t * (smax - smin)


def fmt_p(p) -> str:
    p = float(p)
    return "< 1e-300" if p == 0 else f"{p:.2e}"


def draw_panel(ax, xplot, y, sizes, res, *, log_x, tick_fs, rho_fs, dot_alpha,
               line_lw, ylim=None, show_stats=True):
    log_y = res.get("log_y", False)
    if log_y:                      # zeros have no logarithm and are not plotted
        keep = y > 0
        xplot, y, sizes = xplot[keep], y[keep], sizes[keep]

    order = np.argsort(sizes)
    ax.scatter(xplot[order], y[order], s=sizes[order], c=MARK, alpha=dot_alpha,
               linewidths=0, rasterized=True, zorder=2)

    g = res["xgrid_plot"]
    ax.fill_between(g, res["lo_plot"], res["hi_plot"], color=LINFIT, alpha=0.22,
                    linewidth=0, zorder=3)
    ax.plot(g, res["fit_plot"], color="white", lw=line_lw + 1.1, zorder=4,
            solid_capstyle="round")
    ax.plot(g, res["fit_plot"], color=LINFIT, lw=line_lw, zorder=5,
            solid_capstyle="round")

    if log_x:
        ax.set_xscale("log")
        ax.xaxis.set_major_locator(LogLocator(base=10, numticks=5))
        ax.set_xlim(left=np.nanmin(xplot) * 0.65, right=np.nanmax(xplot) * 1.55)
    else:
        ax.xaxis.set_major_locator(MaxNLocator(nbins=5, prune="both"))
    if log_y:
        ax.set_yscale("log")
        ax.yaxis.set_major_locator(LogLocator(base=10, numticks=5))
        ax.set_ylim(bottom=np.nanmin(y) * 0.65, top=np.nanmax(y) * 1.55)
    else:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5, prune="both"))

    ax.grid(True, color=GRID, lw=0.4, alpha=0.9, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
        ax.spines[side].set_linewidth(0.6)
    ax.tick_params(axis="both", labelsize=tick_fs, length=2.4, width=0.5,
                   pad=2.0, colors=MUTED)
    for lbl in ax.get_xticklabels() + ax.get_yticklabels():
        lbl.set_color(INK)

    if ylim is not None:
        n_above = int((y > ylim).sum())
        ax.set_ylim(bottom=-0.02 * ylim, top=ylim)
        ax.text(0.975, 0.028,
                f"y axis truncated at {ylim:,.0f}; {n_above:,} blocks above,\n"
                f"all retained in the fit",
                transform=ax.transAxes, ha="right", va="bottom",
                fontsize=rho_fs - 0.4, color=MUTED, linespacing=1.3, zorder=9,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none",
                          alpha=0.82))

    if show_stats:
        hc3 = res["lin_hc3"]
        b = float(hc3.params[1])
        lo_ci, hi_ci = hc3.conf_int(alpha=0.05)[1]
        dropped = res.get("n_dropped", 0)
        txt = (f"slope = {b:.4g}\n"
               f"95% CI {lo_ci:.4g} to {hi_ci:.4g}\n"
               f"p = {fmt_p(hc3.pvalues[1])}  (HC3)\n"
               f"$R^2$ = {res['lin'].rsquared:.4f}   n = {res['n']:,}")
        if dropped:
            txt += f"\n{dropped:,} zeros dropped"
        ax.text(0.975, 0.955, txt, transform=ax.transAxes, ha="right", va="top",
                fontsize=rho_fs, color=INK, linespacing=1.35,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none",
                          alpha=0.8), zorder=8)


def legend_handles(smin, smax, df):
    lg = np.log10(np.clip(df[SIZE_COL].to_numpy(float), 1, None))
    lo, hi = np.nanpercentile(lg, 0.5), np.nanpercentile(lg, 99.5)
    handles, labels = [], []
    for L, lab in ((1e4, "10 kb"), (1e5, "100 kb"), (1e6, "1 Mb")):
        t = np.clip((np.log10(L) - lo) / (hi - lo), 0, 1)
        area = smin + t * (smax - smin)
        handles.append(Line2D([], [], marker="o", linestyle="none", color=MARK,
                              markersize=np.sqrt(area), alpha=0.65,
                              markeredgewidth=0))
        labels.append(lab)
    handles.append(Line2D([], [], color=LINFIT, lw=1.6))
    labels.append("OLS fit, 95% CI on the mean")
    return handles, labels


def rcparams(base_fs):
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": base_fs,
        "axes.linewidth": 0.6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def make_figures(df, ycol, resA, resB, resC, outdir: Path):
    x = df[XCOL].to_numpy(float)
    y = df[ycol].to_numpy(float)

    # A zoomed companion row earns its place only when the extreme tail is what
    # compresses the bulk of the data. For a tightly ranged response it would
    # duplicate the panel above it, so it is skipped.
    p99 = float(np.quantile(y, 0.99))
    use_zoom = (y.max() / p99) > ZOOM_RATIO if p99 > 0 else False
    zoom = p99 if use_zoom else None
    RAW, LOGX = "Cluster count", "Cluster count (log scale)"

    # --- combined figure, 180 mm wide -------------------------------------- #
    # Top row: the three model variants at full range. Bottom row: the zoomed
    # companions, which only variants A and B need; the log-log panel already
    # shows its whole range, so its cell is left empty.
    rcparams(7)
    smin, smax = 0.5, 16.0
    sizes = marker_sizes(df[SIZE_COL].to_numpy(), smin, smax)
    nrow = 2 if use_zoom else 1
    fig, axes = plt.subplots(
        nrow, 3, figsize=(180 / 25.4, (152 if use_zoom else 86) / 25.4),
        squeeze=False)
    fig.subplots_adjust(left=0.065, right=0.988, top=0.945,
                        bottom=0.145 if use_zoom else 0.265,
                        wspace=0.30, hspace=0.34)
    cells = [(axes[0][0], resA, RAW, False, "a", None),
             (axes[0][1], resB, LOGX, True, "b", None),
             (axes[0][2], resC, LOGX, True, "c", None)]
    if use_zoom:
        cells += [(axes[1][0], resA, RAW, False, "d", zoom),
                  (axes[1][1], resB, LOGX, True, "e", zoom)]
        axes[1][2].axis("off")
    for ax, res, xlab, logx, letter, ylim in cells:
        draw_panel(ax, x, y, sizes, res, log_x=logx, tick_fs=6, rho_fs=6,
                   dot_alpha=0.06, line_lw=1.4, ylim=ylim)
        ax.set_xlabel(xlab, fontsize=7.5, fontweight="bold", color=INK,
                      labelpad=3)
        ylab = f"{ycol} (log scale)" if res.get("log_y") else ycol
        ax.set_ylabel(ylab, fontsize=7.5, fontweight="bold", color=INK,
                      labelpad=3)
        ax.set_title(letter, loc="left", fontsize=8.5, fontweight="bold",
                     color=INK, pad=3)

    tail = ("; a to c show the full y range, d and e rescale a and b"
            if use_zoom else "")
    handles, labels = legend_handles(smin, smax, df)
    leg = fig.legend(handles, labels, loc="lower center", ncol=len(labels),
                     frameon=False, fontsize=6.2, handletextpad=0.55,
                     columnspacing=1.5, borderaxespad=0.0,
                     bbox_to_anchor=(0.5, 0.008),
                     title=f"block length encoded as dot area; OLS fit{tail}",
                     title_fontsize=6.2)
    leg.get_title().set_color(INK)
    for ext in FORMATS:
        fig.savefig(outdir / f"fig_linear_{ycol}.{ext}", dpi=600,
                    facecolor="white")
    plt.close(fig)

    # --- standalone single-column panels ----------------------------------- #
    rcparams(7)
    smin, smax = 0.5, 15.0
    sizes = marker_sizes(df[SIZE_COL].to_numpy(), smin, smax)
    specs = [(resA, RAW, False, "a", None, "a_raw_x"),
             (resB, LOGX, True, "b", None, "b_log10_x"),
             (resC, LOGX, True, "c", None, "c_loglog")]
    if use_zoom:
        specs += [(resA, RAW, False, "d", zoom, "d_raw_x_zoom"),
                  (resB, LOGX, True, "e", zoom, "e_log10_x_zoom")]
    for res, xlab, logx, letter, ylim, tag in specs:
        fig, ax = plt.subplots(figsize=(89 / 25.4, 84 / 25.4))
        fig.subplots_adjust(left=0.185, right=0.965, top=0.935, bottom=0.315)
        draw_panel(ax, x, y, sizes, res, log_x=logx, tick_fs=6, rho_fs=6,
                   dot_alpha=0.07, line_lw=1.4, ylim=ylim)
        ax.set_xlabel(xlab, fontsize=7.5, fontweight="bold", color=INK,
                      labelpad=3)
        ylab = f"{ycol} (log scale)" if res.get("log_y") else ycol
        ax.set_ylabel(ylab, fontsize=7.5, fontweight="bold", color=INK,
                      labelpad=3)
        ax.set_title(letter, loc="left", fontsize=8, fontweight="bold",
                     color=INK, pad=3)
        handles, labels = legend_handles(smin, smax, df)
        leg = fig.legend(handles, labels, loc="lower center", ncol=3,
                         frameon=False, fontsize=6, handletextpad=0.5,
                         columnspacing=1.2, borderaxespad=0.0,
                         bbox_to_anchor=(0.5, 0.012),
                         title="block length encoded as dot area",
                         title_fontsize=6)
        leg.get_title().set_color(INK)
        for ext in FORMATS:
            fig.savefig(outdir / f"panel_linear_{ycol}_{tag}.{ext}", dpi=600,
                        facecolor="white")
        plt.close(fig)


# --------------------------------------------------------------------------- #
# Report
# --------------------------------------------------------------------------- #

def model_table(w, model, hc3, names):
    w(f"    {'term':<34}{'estimate':>14}{'SE':>12}{'SE (HC3)':>12}"
      f"{'t (HC3)':>10}{'p (HC3)':>12}")
    ci = hc3.conf_int(alpha=0.05)
    for i, nm in enumerate(names):
        w(f"    {nm:<34}{model.params[i]:>14.6g}{model.bse[i]:>12.4g}"
          f"{hc3.bse[i]:>12.4g}{hc3.tvalues[i]:>10.3f}"
          f"{fmt_p(hc3.pvalues[i]):>12}")
        w(f"    {'':<34}{'95% CI (HC3)':>14}  [{ci[i][0]:.6g}, {ci[i][1]:.6g}]")
    w(f"    R2 = {model.rsquared:.6f}   adj R2 = {model.rsquared_adj:.6f}"
      f"   RMSE = {np.sqrt(model.mse_resid):.4f}   df_resid = {int(model.df_resid):,}")
    w("")


def write_report(path: Path, indir: Path, df, fits):
    """fits: list of (ycol, resA, resB), in the order they were modelled."""
    L = []
    w = L.append
    x = df[XCOL].to_numpy(float)

    w("=" * 78)
    w("CLUSTER COUNT vs AVI BLOCK SCORES")
    w("ordinary least squares")
    w("=" * 78)
    w(f"generated        : {_dt.datetime.now().isoformat(timespec='seconds')}")
    w(f"input directory  : {indir}")
    w(f"x                : {XCOL}")
    w(f"y                : {', '.join(y for y, *_ in fits)}")
    w(f"blocks analysed  : {len(df):,}")
    w("")
    w("Variant A regresses raw y on raw x; variant B regresses raw y on")
    w("log10(x); variant C is the log-log model, log10(y) on log10(x).")
    w("")

    w("-" * 78)
    w("1. DISTRIBUTIONS")
    w("-" * 78)
    w(f"  {'':<20}{'min':>12}{'p10':>12}{'median':>12}{'p90':>12}{'max':>14}")
    rows = [(XCOL, x)] + [(yc, df[yc].to_numpy(float)) for yc, *_ in fits]
    rows.append((SIZE_COL, df[SIZE_COL].to_numpy(float)))
    for nm, v in rows:
        q = np.quantile(v, [0, .10, .50, .90, 1])
        w(f"  {nm:<20}{q[0]:>12,.4g}{q[1]:>12,.4g}{q[2]:>12,.4g}"
          f"{q[3]:>12,.4g}{q[4]:>14,.6g}")
    w("")
    for yc, *_ in fits:
        y = df[yc].to_numpy(float)
        w(f"  {yc}")
        w(f"    blocks == 0                : {(y == 0).sum():,} "
          f"({100*(y == 0).mean():.2f}%), retained because y is untransformed")
        w(f"    Spearman rho (x, y)        : {stats.spearmanr(x, y)[0]:.4f}")
        w(f"    Pearson r (raw x, raw y)   : {stats.pearsonr(x, y)[0]:.4f}")
        w(f"    Pearson r (log10 x, raw y) : "
          f"{stats.pearsonr(np.log10(x), y)[0]:.4f}")
    w("")

    sections = []
    for yc, rA, rB, rC in fits:
        sections.append((yc, rA, f"VARIANT A - raw {yc} on raw x", XCOL))
        sections.append((yc, rB, f"VARIANT B - raw {yc} on log10 x",
                         f"log10({XCOL})"))
        sections.append((yc, rC, f"VARIANT C - log10 {yc} on log10 x",
                         f"log10({XCOL})"))

    n = 1
    for ycol, res, title, xdesc in sections:
        n += 1
        w("-" * 78)
        w(f"{n}. {title}")
        w("-" * 78)
        w(f"  predictor : {xdesc}")
        if res.get("log_y"):
            w(f"  response  : log10({ycol})")
        else:
            w(f"  response  : {ycol} (untransformed)")
        w(f"  n         : {res['n']:,}")
        if res.get("n_dropped"):
            w(f"  dropped   : {res['n_dropped']:,} blocks with {ycol} == 0, "
              "which have no logarithm")
        if res.get("log_y"):
            w("  the slope is an elasticity: the proportional change in y for a")
            w("  proportional change in x. A slope of 1 means y doubles when x")
            w("  doubles.")
        w("")
        ydesc = f"log10({ycol})" if res.get("log_y") else ycol
        model_table(w, res["lin"], res["lin_hc3"], ["intercept", xdesc])

    w("-" * 78)
    w(f"{n + 1}. NOTES ON INTERPRETATION")
    w("-" * 78)
    w("  With n near 39,000 every coefficient here is significant at any")
    w("  conventional threshold. Read the effect sizes and R2, not the")
    w("  p-values.")
    w("")
    w("  HC3 standard errors are reported alongside the classical ones because")
    w("  the residual variance rises steeply with x; the classical SEs are")
    w("  optimistic. Quote the HC3 column.")
    w("")
    w("  Variant A puts a small number of very high cluster-count blocks at")
    w("  extreme leverage: the top 1% of x reaches "
      f"{np.quantile(x, 0.99):,.0f} against a median of {np.median(x):,.0f}.")
    w("  Variant B is the better-conditioned fit for that reason.")
    w("")
    w("  The 95% bands are on the conditional mean, not prediction intervals.")
    w("  At this sample size they are narrow enough to be nearly invisible on")
    w("  the figure; that reflects precision about the mean, not about any")
    w("  individual block.")
    w("")
    w("=" * 78)
    w("END OF REPORT")
    w("=" * 78)
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


def write_fitted(path: Path, fits):
    rows = []
    specs = []
    for ycol, rA, rB, rC in fits:
        specs += [(ycol, rA, "A_raw_x"), (ycol, rB, "B_log10_x"),
                  (ycol, rC, "C_loglog")]
    for ycol, res, variant in specs:
        for i in range(len(res["xgrid"])):
            rows.append({
                "y": ycol,
                "variant": variant,
                "y_transform": "log10" if res.get("log_y") else "raw",
                "x_model": res["xgrid"][i],
                "x_plot": res["xgrid_plot"][i],
                "fit_model_space": res["lin_fit"][i],
                "ci_lo_model_space": res["lin_lo"][i],
                "ci_hi_model_space": res["lin_hi"][i],
                "fit_y_units": res["fit_plot"][i],
                "ci_lo_y_units": res["lo_plot"][i],
                "ci_hi_y_units": res["hi_plot"][i],
            })
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, float_format="%.6g")


# --------------------------------------------------------------------------- #

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", default=".", type=Path)
    ap.add_argument("--outdir", default=".", type=Path)
    args = ap.parse_args(argv)

    indir, outdir = args.indir.resolve(), args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    df = load_and_join(indir)
    missing = [c for c in Y_TARGETS if c not in df.columns]
    if missing:
        raise SystemExit(f"missing response column(s) in {AVI_FILE}: {missing}")

    keep = [XCOL, SIZE_COL, "_key"] + Y_TARGETS
    df = df[keep].replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=[XCOL, SIZE_COL] + Y_TARGETS)
    df = df[df[XCOL] > 0]           # log10 needs positive x; min is 1 anyway

    x = df[XCOL].to_numpy(float)
    lx = np.log10(x)
    gridA = np.linspace(x.min(), x.max(), 400)
    gridB = np.linspace(lx.min(), lx.max(), 400)

    fits = []
    for ycol in Y_TARGETS:
        y = df[ycol].to_numpy(float)

        resA = fit_variant(x, y, gridA)          # raw x, grid linear in x
        resA["xgrid_plot"] = gridA

        resB = fit_variant(lx, y, gridB)         # log10 x, plotted on a log axis
        resB["xgrid_plot"] = 10 ** gridB

        # Variant C, log-log. Zeros have no logarithm, so those blocks are
        # dropped from this fit alone; the count is carried into the report and
        # printed on the panel.
        pos = y > 0
        lxc, lyc = lx[pos], np.log10(y[pos])
        gridC = np.linspace(lxc.min(), lxc.max(), 400)
        resC = fit_variant(lxc, lyc, gridC, log_y=True)
        resC["xgrid_plot"] = 10 ** gridC
        resC["n_dropped"] = int((~pos).sum())

        make_figures(df, ycol, resA, resB, resC, outdir)
        fits.append((ycol, resA, resB, resC))

    write_report(outdir / "regression_report.txt", indir, df, fits)
    write_fitted(outdir / "fitted_values.tsv", fits)

    print(f"n = {len(df):,}")
    for ycol, rA, rB, rC in fits:
        print(f"{ycol}:  A R2 = {rA['lin'].rsquared:.5f}   "
              f"B R2 = {rB['lin'].rsquared:.5f}   "
              f"C R2 = {rC['lin'].rsquared:.5f} (elasticity "
              f"{rC['lin_hc3'].params[1]:.4f}, n dropped {rC['n_dropped']})")
    print(f"wrote outputs to {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

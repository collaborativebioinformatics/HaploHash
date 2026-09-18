"""Test block AVI enrichment in GTEx v8 using distinct tested SNVs as background."""

from __future__ import annotations

import argparse
import os
import re
import shutil
import tarfile
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE / "gtex"
DATA, WORK, RESULTS = ROOT / "data", ROOT / "work", ROOT / "results"
OUTCOMES = {"eqtl": "is_eqtl", "lead": "is_lead", "finemapped": "is_finemapped"}
MAF_BINS = [0.02, 0.05, 0.1, 0.2, 0.3, 0.4]
DISTANCE_BINS = [1000, 5000, 10000, 25000, 50000, 100000, 250000, 500000]
GENE_BINS = [2, 4, 8, 16, 32]


def sql_string(value):
    return "'" + str(value).replace("'", "''") + "'"


def connection(name):
    temp = WORK / "tmp" / name
    temp.mkdir(parents=True, exist_ok=True)
    return duckdb.connect(config={
        "threads": int(os.environ.get("SLURM_CPUS_PER_TASK", 4)),
        "memory_limit": "18GB", "temp_directory": str(temp),
        "preserve_insertion_order": False,
    })


def tissue_names():
    return pd.read_csv(ROOT / "tissues.tsv", sep="\t").tissue.tolist()


def canonical_tissue(name, names):
    normalize = lambda s: re.sub(r"[^a-z0-9]", "", s.lower())
    key = normalize(name)
    matches = [s for s in names if normalize(s) == key]
    if not matches:
        matches = [s for s in names if normalize(s).startswith(key)]
    if len(matches) != 1:
        raise ValueError(f"ambiguous fine-mapping tissue {name!r}: {matches}")
    return matches[0]


def prepare():
    WORK.mkdir(parents=True, exist_ok=True)
    archive_path = DATA / "GTEx_v8_finemapping_DAPG.tar"
    source = WORK / "dapg.txt.gz"
    if not source.exists():
        with tarfile.open(archive_path) as archive:
            member = archive.getmember("GTEx_v8_finemapping_DAPG/GTEx_v8_finemapping_DAPG.txt.gz")
            with archive.extractfile(member) as src, source.open("wb") as dst:
                shutil.copyfileobj(src, dst)
    con = connection("prepare")
    fine = con.execute("""
        SELECT tissue_id AS tissue, regexp_replace(variant_id, '_b38$', '') AS variant,
               max(variant_pip) AS pip
        FROM read_csv(?, delim='\t', header=true)
        WHERE variant_pip >= 0.5
        GROUP BY tissue, variant
    """, [str(source)]).df()
    names = tissue_names()
    mapping = {t: canonical_tissue(t, names) for t in fine.tissue.unique()}
    fine["tissue"] = fine.tissue.map(mapping)
    (WORK / "finemapped").mkdir(exist_ok=True)
    for tissue in names:
        fine.loc[fine.tissue == tissue, ["variant", "pip"]].to_parquet(
            WORK / "finemapped" / f"{tissue}.parquet", index=False, compression="zstd")
    con.close()
    print(f"Prepared PIP ≥ 0.5 variants for {len(names)} tissues", flush=True)


def load_blocks():
    blocks = pd.read_csv(HERE / "block_scores.tsv", sep="\t")
    blocks = blocks[blocks.chrom.isin([f"chr{i}" for i in range(1, 23)]) & (blocks.n_scored > 0)].copy()
    blocks["avi_fraction"] = blocks.avi_top1_count / blocks.n_scored
    blocks["avi_quintile"] = pd.qcut(blocks.avi_fraction, 5, labels=False) + 1
    blocks["region"] = blocks.chrom + ":" + ((blocks.start + blocks.end) // 10_000_000).astype(str)
    blocks = blocks.sort_values(["chrom", "start"]).reset_index(drop=True)
    blocks["block_index"] = np.arange(len(blocks))
    return blocks


def assign_blocks(variants, blocks):
    assigned = np.full(len(variants), -1, dtype=np.int32)
    for chrom, indices in variants.groupby("chrom", sort=False).indices.items():
        intervals = blocks[blocks.chrom == chrom]
        if intervals.empty:
            continue
        pos = variants.position.to_numpy()[indices] - 1
        starts, ends = intervals.start.to_numpy(), intervals.end.to_numpy()
        j = np.searchsorted(starts, pos, side="right") - 1
        valid = (j >= 0) & (pos < ends[np.maximum(j, 0)])
        assigned[indices[valid]] = intervals.block_index.to_numpy()[j[valid]]
    return assigned


def load_genes(tissue):
    genes = pd.read_csv(DATA / "egenes" / f"{tissue}.v8.egenes.txt.gz", sep="\t")
    genes["gene_id"] = genes.gene_id.str.split(".").str[0]
    if genes.gene_id.duplicated().any():
        raise ValueError("gene IDs are not unique")
    genes["tss"] = np.where(genes.strand == "+", genes.gene_start, genes.gene_end)
    genes["lead_variant"] = genes.variant_id.str.replace("_b38", "", regex=False)
    genes["egene"] = genes.qval <= 0.05
    return genes


def validate_hits(tissue, variants):
    """Use official calls: rounded nominal thresholds can mislabel boundary cases."""
    member = f"GTEx_Analysis_v8_eQTL/{tissue}.v8.signif_variant_gene_pairs.txt.gz"
    with tarfile.open(DATA / "GTEx_Analysis_v8_eQTL.tar") as archive:
        with archive.extractfile(member) as stream:
            ids = pd.read_csv(stream, compression="gzip", sep="\t", usecols=["variant_id"]).variant_id.drop_duplicates()
    ids = ids[ids.str.match(r"chr(?:[1-9]|1[0-9]|2[0-2])_\d+_[ACGT]_[ACGT]_b38$")]
    official = set(ids.str.replace("_b38", "", regex=False))
    observed = set(variants.loc[variants.is_eqtl, "variant"])
    missing = official.difference(variants.variant)
    if missing:
        raise ValueError(f"{tissue}: {len(missing)} official eQTL SNVs absent from tested background")
    if official != observed:
        print(f"{tissue}: using official calls for {len(official ^ observed)} threshold-label discrepancies", flush=True)
    variants["is_eqtl"] = variants.variant.isin(official)
    variants["is_lead"] &= variants.is_eqtl


def collapse_variants(tissue, variants_path):
    """Parse once, then aggregate chromosomes separately to bound memory use."""
    pairs_path = WORK / f"{tissue}.pairs.parquet"
    con = connection(tissue)
    if not pairs_path.exists():
        print(f"{tissue}: parsing all tested variant–gene pairs", flush=True)
        genes = load_genes(tissue)
        con.register("genes", genes)
        source = DATA / "allpairs" / f"{tissue}.tsv.gz"
        con.execute(f"""
            COPY (
                SELECT s.variant, 'chr' || s.chromosome AS chrom,
                       CAST(s.position AS BIGINT) AS position,
                       CAST(s.maf AS DOUBLE) AS maf, s.gene_id,
                       abs(CAST(s.position AS BIGINT) - g.tss) AS min_tss_distance,
                       g.egene AND CAST(s.pvalue AS DOUBLE) <= g.pval_nominal_threshold AS is_eqtl,
                       g.egene AND s.variant = g.lead_variant AS is_lead,
                       g.gene_id IS NULL AS missing_gene,
                       CAST(s.an AS BIGINT) AS an
                FROM read_csv({sql_string(source)}, delim='\t', header=true,
                              all_varchar=true, nullstr=['NA', 'NaN']) s
                LEFT JOIN genes g ON s.gene_id = g.gene_id
                WHERE s.type = 'SNP' AND length(s.ref) = 1 AND length(s.alt) = 1
                      AND s.chromosome IN ({','.join(sql_string(i) for i in range(1, 23))})
            ) TO {sql_string(str(pairs_path) + '.part')} (FORMAT PARQUET, COMPRESSION ZSTD)
        """)
        Path(str(pairs_path) + ".part").rename(pairs_path)
    pieces = WORK / f"{tissue}.chromosomes"
    pieces.mkdir(exist_ok=True)
    for chrom in range(1, 23):
        target = pieces / f"chr{chrom}.parquet"
        if target.exists():
            continue
        # Collapse duplicated rsIDs and all gene tests to one row per alternate allele.
        con.execute(f"""
            COPY (
                SELECT variant, chrom, position, min(maf) AS maf,
                       count(DISTINCT gene_id) AS n_tested_genes,
                       min(min_tss_distance) AS min_tss_distance,
                       bool_or(is_eqtl) AS is_eqtl, bool_or(is_lead) AS is_lead,
                       bool_or(missing_gene) AS missing_gene, max(an) AS an
                FROM read_parquet({sql_string(pairs_path)})
                WHERE chrom = 'chr{chrom}'
                GROUP BY variant, chrom, position
            ) TO {sql_string(str(target) + '.part')} (FORMAT PARQUET, COMPRESSION ZSTD)
        """)
        Path(str(target) + ".part").rename(target)
        print(f"{tissue}: collapsed chr{chrom}", flush=True)
    con.execute(f"COPY (SELECT * FROM read_parquet({sql_string(pieces / '*.parquet')})) "
                f"TO {sql_string(str(variants_path) + '.part')} (FORMAT PARQUET, COMPRESSION ZSTD)")
    con.close()
    Path(str(variants_path) + ".part").rename(variants_path)
    pairs_path.unlink()
    shutil.rmtree(pieces)


def process_tissue(tissue):
    WORK.mkdir(parents=True, exist_ok=True)
    variants_path = WORK / f"{tissue}.variants.parquet"
    if not variants_path.exists():
        collapse_variants(tissue, variants_path)
    variants = pd.read_parquet(variants_path)
    if variants.missing_gene.any():
        raise ValueError(f"{tissue}: variants have tests for genes absent from official eGene table")
    if variants[["maf", "min_tss_distance", "is_eqtl", "is_lead"]].isna().any().any():
        raise ValueError(f"{tissue}: missing comparison data")
    if not variants.maf.between(0, 0.5).all():
        raise ValueError(f"{tissue}: invalid MAF")
    validate_hits(tissue, variants)
    print(f"{tissue}: {len(variants):,} tested SNVs; official GTEx calls covered", flush=True)
    fine = pd.read_parquet(WORK / "finemapped" / f"{tissue}.parquet")
    variants["is_finemapped"] = variants.variant.isin(fine.variant) & variants.is_eqtl
    blocks = load_blocks()
    variants["block_index"] = assign_blocks(variants, blocks)
    n_eligible = len(variants)
    sample_size = int(variants.an.max() // 2)
    variants = variants[variants.block_index >= 0].copy()
    by_stratum = stratify(variants)
    by_stratum.to_parquet(WORK / f"{tissue}.strata.parquet", index=False, compression="zstd")
    print(f"{tissue}: {len(variants):,} SNVs in blocks; calculating regional intervals", flush=True)
    block_counts, quintiles, summary = analyze_strata(by_stratum, blocks, tissue)
    summary["sample_size"] = sample_size
    summary["tested_autosomal_snvs"] = n_eligible
    block_counts.to_parquet(WORK / f"{tissue}.blocks.parquet", index=False, compression="zstd")
    quintiles.to_csv(WORK / f"{tissue}.quintiles.tsv", sep="\t", index=False)
    summary.to_csv(WORK / f"{tissue}.summary.tsv", sep="\t", index=False)
    print(summary[["tissue", "outcome", "raw_high_low_ratio", "adjusted_high_low_ratio", "ci_low", "ci_high"]].to_string(index=False), flush=True)


def stratify(variants, maf_bins=MAF_BINS, distance_bins=DISTANCE_BINS, gene_bins=GENE_BINS):
    maf = np.searchsorted(maf_bins, variants.maf, side="right")
    distance = np.searchsorted(distance_bins, variants.min_tss_distance, side="right")
    genes = np.searchsorted(gene_bins, variants.n_tested_genes, side="right")
    strata = (maf * (len(distance_bins) + 1) + distance) * (len(gene_bins) + 1) + genes
    return variants.groupby([variants.block_index, pd.Series(strata, index=variants.index, name="stratum")], sort=False).agg(
        tested=("variant", "size"), eqtl=("is_eqtl", "sum"),
        lead=("is_lead", "sum"), finemapped=("is_finemapped", "sum"),
    ).reset_index()


def analyze_strata(strata, blocks, tissue, n_bootstrap=1000, seed=20260917):
    """Standardize each block to tissue-wide outcome rates within covariate strata."""
    totals = strata.groupby("stratum")[["tested", *OUTCOMES]].sum()
    counts = strata.groupby("block_index")[["tested", *OUTCOMES]].sum()
    counts = blocks.merge(counts, left_on="block_index", right_index=True, how="left")
    counts[["tested", *OUTCOMES]] = counts[["tested", *OUTCOMES]].fillna(0).astype(np.int64)
    rng = np.random.default_rng(seed)
    regions = sorted(counts.region.unique())
    # Resample 5 Mb genomic regions, keeping correlated variants/blocks together.
    weights = rng.multinomial(len(regions), np.full(len(regions), 1 / len(regions)), size=n_bootstrap)
    regional = strata.merge(blocks[["block_index", "region", "avi_quintile"]], on="block_index", validate="many_to_one")
    cells = pd.MultiIndex.from_product([totals.index, range(1, 6)], names=["stratum", "avi_quintile"])

    def resample(column):
        matrix = regional.pivot_table(index="region", columns=["stratum", "avi_quintile"],
                                      values=column, aggfunc="sum", fill_value=0)
        matrix = matrix.reindex(index=regions, columns=cells, fill_value=0).fillna(0)
        return (weights @ matrix.to_numpy(dtype=float)).reshape(n_bootstrap, len(totals), 5)

    boot_tested = resample("tested")
    quintiles, summary = [], []
    for outcome in OUTCOMES:
        expected = strata.tested * strata.stratum.map(totals[outcome] / totals.tested)
        expected = expected.groupby(strata.block_index).sum()
        counts[f"expected_{outcome}"] = counts.block_index.map(expected).fillna(0)
        pooled = counts.groupby("avi_quintile")[["tested", outcome, f"expected_{outcome}"]].sum().reindex(range(1, 6), fill_value=0)
        boot_cases = resample(outcome)
        # Refit each stratum's reference rate within every bootstrap replicate.
        boot_rate = np.divide(boot_cases.sum(axis=2), boot_tested.sum(axis=2),
                              out=np.zeros(boot_cases.shape[:2]), where=boot_tested.sum(axis=2) > 0)
        boot_observed = boot_cases.sum(axis=1)
        boot_expected = (boot_tested * boot_rate[:, :, None]).sum(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            boot_oe = boot_observed / boot_expected
            boot_ratio = boot_oe[:, 4] / boot_oe[:, 0]
            oe = pooled[outcome].to_numpy() / pooled[f"expected_{outcome}"].to_numpy()
        for i, (q, row) in enumerate(pooled.iterrows()):
            values = boot_oe[:, i]
            values = values[np.isfinite(values)]
            low, high = np.quantile(values, [0.025, 0.975]) if len(values) else (np.nan, np.nan)
            quintiles.append(dict(tissue=tissue, outcome=outcome, avi_quintile=q,
                                  tested=int(row.tested), observed=int(row[outcome]),
                                  expected=row[f"expected_{outcome}"], observed_expected=oe[i],
                                  ci_low=low, ci_high=high))
        values = boot_ratio[np.isfinite(boot_ratio)]
        low, high = np.quantile(values, [0.025, 0.975]) if len(values) else (np.nan, np.nan)
        rates = pooled[outcome] / pooled.tested.replace(0, np.nan)
        raw = rates.loc[5] / rates.loc[1] if rates.loc[1] > 0 else np.nan
        ratio = oe[4] / oe[0] if oe[0] > 0 else np.nan
        summary.append(dict(tissue=tissue, outcome=outcome,
                            n_tested=int(counts.tested.sum()), n_observed=int(counts[outcome].sum()),
                            raw_high_low_ratio=raw, adjusted_high_low_ratio=ratio,
                            ci_low=low, ci_high=high))
    counts["tissue"] = tissue
    return counts, pd.DataFrame(quintiles), pd.DataFrame(summary)


def combine():
    names = tissue_names()
    missing = [t for t in names if not (WORK / f"{t}.summary.tsv").exists()]
    if missing:
        raise ValueError(f"unfinished tissues: {', '.join(missing)}")
    RESULTS.mkdir(parents=True, exist_ok=True)
    summary = pd.concat([pd.read_csv(WORK / f"{t}.summary.tsv", sep="\t") for t in names], ignore_index=True)
    quintiles = pd.concat([pd.read_csv(WORK / f"{t}.quintiles.tsv", sep="\t") for t in names], ignore_index=True)
    summary.to_csv(RESULTS / "tissue_summary.tsv", sep="\t", index=False, float_format="%.6g")
    quintiles.to_csv(RESULTS / "avi_quintiles.tsv", sep="\t", index=False, float_format="%.6g")
    columns = ["tissue", "block_id", "chrom", "start", "end", "avi_fraction", "avi_quintile",
               "tested", *OUTCOMES, *[f"expected_{k}" for k in OUTCOMES]]
    all_blocks = pd.concat([pd.read_parquet(WORK / f"{t}.blocks.parquet", columns=columns) for t in names], ignore_index=True)
    all_blocks.to_parquet(RESULTS / "block_eqtl_counts.parquet", compression="zstd", index=False)
    del all_blocks
    blood_sensitivity()
    plot_results(summary, quintiles)
    print(summary[summary.tissue == "Whole_Blood"].to_string(index=False))


def blood_sensitivity():
    tissue = "Whole_Blood"
    blocks = load_blocks()
    baseline = pd.read_csv(WORK / f"{tissue}.summary.tsv", sep="\t")
    baseline["comparison"] = "primary"
    strata = pd.read_parquet(WORK / f"{tissue}.strata.parquet")
    subset = blocks[blocks.chrom != "chr6"]
    _, _, no_chr6 = analyze_strata(strata[strata.block_index.isin(subset.block_index)], subset, tissue)
    no_chr6["comparison"] = "exclude_chr6"
    variants = pd.read_parquet(WORK / f"{tissue}.variants.parquet")
    validate_hits(tissue, variants)
    fine = pd.read_parquet(WORK / "finemapped" / f"{tissue}.parquet")
    variants["is_finemapped"] = variants.variant.isin(fine.variant) & variants.is_eqtl
    variants["block_index"] = assign_blocks(variants, blocks)
    variants = variants[variants.block_index >= 0]
    finer = stratify(variants,
                    maf_bins=[0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45],
                    distance_bins=[500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 250000, 500000],
                    gene_bins=[2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64])
    del variants
    _, _, finer_summary = analyze_strata(finer, blocks, tissue)
    finer_summary["comparison"] = "finer_covariate_bins"
    columns = ["comparison", "outcome", "n_tested", "n_observed", "adjusted_high_low_ratio", "ci_low", "ci_high"]
    result = pd.concat([baseline, no_chr6, finer_summary], ignore_index=True)[columns]
    RESULTS.mkdir(parents=True, exist_ok=True)
    result.to_csv(RESULTS / "blood_sensitivity.tsv", sep="\t", index=False, float_format="%.6g")
    print(result.to_string(index=False), flush=True)


def plot_results(summary, quintiles):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter

    plt.style.use(HERE / "paper.mplstyle")
    labels = {"eqtl": "eQTLs", "lead": "Lead eQTLs", "finemapped": "Fine-mapped (PIP ≥ 0.5)"}
    colors = dict(zip(OUTCOMES, plt.rcParams["axes.prop_cycle"].by_key()["color"]))
    markers = {"eqtl": "o", "lead": "s", "finemapped": "^"}
    blood = quintiles[quintiles.tissue == "Whole_Blood"]
    fig, axes = plt.subplots(1, 2, figsize=(7, 2.7), layout="constrained")
    for offset, outcome in zip([-0.055, 0, 0.055], OUTCOMES):
        d = blood[blood.outcome == outcome].sort_values("avi_quintile")
        axes[0].plot(d.avi_quintile, d.observed / d.tested, marker=markers[outcome],
                     color=colors[outcome], label=labels[outcome])
        axes[1].errorbar(d.avi_quintile + offset, d.observed_expected,
                         yerr=[np.maximum(0, d.observed_expected - d.ci_low), np.maximum(0, d.ci_high - d.observed_expected)],
                         marker=markers[outcome], capsize=1.8, capthick=0.7, elinewidth=0.7,
                         color=colors[outcome])
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Fraction of tested SNVs")
    axes[0].set_title("GTEx v8 whole blood: unadjusted")
    axes[0].legend(loc="lower right", fontsize=6.5)
    axes[1].set_ylabel("Observed / expected (95% CI)")
    axes[1].set_title("Covariate-adjusted")
    axes[1].axhline(1, color="0.4", linestyle="--", lw=0.7)
    for ax in axes:
        ax.set_xticks(range(1, 6))
        ax.set_xlim(0.75, 5.25)
        ax.set_xlabel("AVI quintile (low → high)")
    fig.savefig(RESULTS / "whole_blood.png")
    plt.close(fig)

    order = summary[summary.outcome == "eqtl"].sort_values("adjusted_high_low_ratio").tissue.tolist()
    fig, axes = plt.subplots(1, 3, figsize=(8, 9.5), sharey=True, layout="constrained")
    for ax, outcome in zip(axes, OUTCOMES):
        d = summary[summary.outcome == outcome].set_index("tissue").loc[order]
        for i, (_, row) in enumerate(d.iterrows()):
            is_blood = row.name == "Whole_Blood"
            color = colors["eqtl"] if is_blood else "#333333"
            if is_blood:
                ax.axhspan(i - 0.45, i + 0.45, color=colors["eqtl"], alpha=0.06, lw=0)
            ax.plot([row.ci_low, row.ci_high], [i, i], color=color, lw=0.7)
            ax.scatter(row.adjusted_high_low_ratio, i, color=color, s=12 if is_blood else 7, zorder=3)
        ax.axvline(1, color="0.45", linestyle="--", lw=0.7)
        ax.set_xscale("log")
        low, high = ax.get_xlim()
        ticks = [t for t in [0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 10] if low <= t <= high]
        ax.set_xticks(ticks, [f"{t:g}" for t in ticks])
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_title("GTEx v8: eQTLs" if outcome == "eqtl" else labels[outcome], fontsize=8)
        ax.grid(axis="y", visible=False)
        ax.tick_params(axis="y", length=0)
        ax.margins(y=0.012)
    axes[0].set_yticks(range(len(order)), [t.replace("_", " ") for t in order], fontsize=7)
    axes[0].invert_yaxis()
    for label in axes[0].get_yticklabels():
        if label.get_text() == "Whole Blood":
            label.set_color(colors["eqtl"])
            label.set_fontweight("bold")
    fig.supxlabel("Adjusted enrichment: highest / lowest AVI quintile (95% CI)", fontsize=8)
    fig.savefig(RESULTS / "across_tissues.png")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["prepare", "tissue", "combine", "sensitivity"])
    parser.add_argument("--index", type=int, default=0)
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare()
    elif args.stage == "tissue":
        process_tissue(tissue_names()[args.index])
    elif args.stage == "sensitivity":
        blood_sensitivity()
    else:
        combine()


if __name__ == "__main__":
    main()

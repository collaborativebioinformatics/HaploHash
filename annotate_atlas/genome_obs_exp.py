#!/usr/bin/env python3
"""Genome-wide obs/exp constraint: for every haploblock, how many of its
AVI top1-tier (PHRED>=20) positions are actually observed as variants in
1000G, versus how many such positions exist.

Generalizes build_variant_parquet.py's per-variant join (12 chr22 pilot
blocks only) to all 39,074 blocks. Streams and aggregates directly to
per-block counts instead of materializing every variant, since genome-wide
this is tens of millions of rows -- a full parquet isn't needed for this
ratio.

obs_top1_count / exp_top1_count (exp_top1_count = avi_top1_count column in
block_scores.tsv) is the constraint ratio: low = most potentially-damaging
positions never show up as real variation (selection is suppressing them).

Needs one phased VCF per chromosome (only chr22's is on the cluster right
now -- see data/1000G/). Contigs are assumed bare ("22" not "chr22"), matching
the existing 1000G download; --contig-style lets that be overridden.

Usage:
  python3 genome_obs_exp.py blocks.bed VCF_DIR AVI_PATH obs_exp_by_block.tsv \
      --tabix path/to/tabix --bcftools bcftools
VCF_DIR must contain one bgzipped+tabixed VCF per chromosome, named so that
glob "*.vcf.gz" resolves to exactly one file per contig (matched via
`bcftools query -r <contig>` against each file's own contig list).
"""

from __future__ import annotations

import argparse
import bisect
import glob
import os
import subprocess
import sys
from dataclasses import dataclass, field

sys.path.insert(0, os.path.dirname(__file__))
from score_blocks import (
    canonical_chrom,
    detect_contig_style,
    parse_column_spec,
    read_tabix_contigs,
    read_tabix_header,
    tabix_chrom,
)
from build_variant_parquet import (
    AVI_ALT_ALIASES,
    AVI_PHRED_ALIASES,
    AVI_POSITION_ALIASES,
)


@dataclass
class Block:
    chrom: str
    start: int
    end: int
    block_id: str
    observed_total: int = 0
    observed_top1: int = 0


def load_blocks(bed_path: str) -> dict[str, list[Block]]:
    by_chrom: dict[str, list[Block]] = {}
    with open(bed_path) as fh:
        for line in fh:
            chrom, start, end, block_id = line.rstrip("\n").split("\t")
            by_chrom.setdefault(chrom, []).append(Block(chrom, int(start), int(end), block_id))
    for blocks in by_chrom.values():
        blocks.sort(key=lambda b: b.start)
    return by_chrom


def vcf_contigs(bcftools: str, vcf_path: str) -> set[str]:
    out = subprocess.run([bcftools, "index", "-s", vcf_path], capture_output=True, text=True, check=True)
    return {line.split("\t")[0] for line in out.stdout.splitlines() if line.strip()}


def find_vcf_for_chrom(bcftools: str, vcf_dir: str, chrom_bare: str) -> str | None:
    for path in glob.glob(os.path.join(vcf_dir, "*.vcf.gz")):
        if chrom_bare in vcf_contigs(bcftools, path):
            return path
    return None


def build_avi_index_for_chrom(tabix: str, avi_path: str, chrom_avi: str,
                               phred_column: int, pos_column: int, alt_column: int,
                               ref_column: int, block: Block) -> dict[tuple[int, str], float]:
    region = f"{chrom_avi}:{block.start + 1}-{block.end}"
    process = subprocess.Popen([tabix, avi_path, region], stdout=subprocess.PIPE, text=True)
    index: dict[tuple[int, str], float] = {}
    assert process.stdout is not None
    for line in process.stdout:
        fields = line.rstrip("\n").split("\t")
        pos = int(fields[pos_column])
        alt = fields[alt_column]
        phred = float(fields[phred_column])
        index[(pos, alt)] = phred
    process.wait()
    return index


def main(argv: list[str]) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("blocks_bed")
    p.add_argument("vcf_dir")
    p.add_argument("avi_path")
    p.add_argument("output")
    p.add_argument("--tabix", default="tabix")
    p.add_argument("--bcftools", default="bcftools")
    p.add_argument("--avi-phred-column", default=None)
    p.add_argument("--avi-position-column", default=None)
    p.add_argument("--avi-alt-column", default=None)
    p.add_argument("--contig-style", choices=("auto", "chr", "bare"), default="auto")
    p.add_argument("--chrom", default=None,
                    help="process only this chromosome (e.g. chr7); for one-process-per-chromosome "
                         "parallelism, matching how avi_score_blocks.sbatch parallelized the static "
                         "genome-wide AVI scoring run")
    args = p.parse_args(argv[1:])

    by_chrom = load_blocks(args.blocks_bed)
    if args.chrom is not None:
        by_chrom = {args.chrom: by_chrom[args.chrom]} if args.chrom in by_chrom else {}

    avi_header = read_tabix_header(args.tabix, args.avi_path)
    contig_style = args.contig_style
    if contig_style == "auto":
        contig_style = detect_contig_style(read_tabix_contigs(args.tabix, args.avi_path))
    pos_col = parse_column_spec(args.avi_position_column, avi_header, AVI_POSITION_ALIASES,
                                 "avi-position", fallback=1)
    alt_col = parse_column_spec(args.avi_alt_column, avi_header, AVI_ALT_ALIASES, "avi-alt")
    phred_col = parse_column_spec(args.avi_phred_column, avi_header, AVI_PHRED_ALIASES, "avi-phred")

    for chrom, blocks in sorted(by_chrom.items()):
        chrom_bare = chrom.removeprefix("chr")
        vcf_path = find_vcf_for_chrom(args.bcftools, args.vcf_dir, chrom_bare)
        if vcf_path is None:
            print(f"[{chrom}] no VCF found in {args.vcf_dir}, skipping", file=sys.stderr)
            continue

        chrom_avi = tabix_chrom(chrom, contig_style)
        starts = [b.start for b in blocks]

        regions_arg = f"{chrom_bare}"
        query = subprocess.Popen(
            [args.bcftools, "query", "-r", regions_arg, "-f", "%POS\t%REF\t%ALT\t%INFO/VT\n", vcf_path],
            stdout=subprocess.PIPE, text=True,
        )
        assert query.stdout is not None

        cursor = 0
        n_lines = 0
        avi_cache: dict[tuple[int, str], float] | None = None
        cache_block_idx = -1
        for line in query.stdout:
            n_lines += 1
            fields = line.rstrip("\n").split("\t")
            pos, ref, alt, vt = fields[0], fields[1], fields[2], fields[3]
            if vt != "SNP" or len(ref) != 1 or len(alt) != 1:
                continue
            pos_int = int(pos)
            pos0 = pos_int - 1
            while cursor < len(blocks) and blocks[cursor].end <= pos0:
                cursor += 1
            if cursor >= len(blocks) or not (blocks[cursor].start <= pos0 < blocks[cursor].end):
                continue
            block = blocks[cursor]

            if cursor != cache_block_idx:
                avi_cache = build_avi_index_for_chrom(
                    args.tabix, args.avi_path, chrom_avi, phred_col, pos_col, alt_col, None, block
                )
                cache_block_idx = cursor

            phred = avi_cache.get((pos_int, alt)) if avi_cache else None
            if phred is None:
                continue
            block.observed_total += 1
            if phred >= 20:
                block.observed_top1 += 1
        query.wait()
        print(f"[{chrom}] {n_lines} VCF rows scanned, "
              f"{sum(b.observed_total for b in blocks)} matched to AVI", file=sys.stderr)

    with open(args.output, "w") as fh:
        fh.write("block_id\tchrom\tstart\tend\tobserved_total\tobserved_top1\n")
        for blocks in by_chrom.values():
            for b in blocks:
                fh.write(f"{b.block_id}\t{b.chrom}\t{b.start}\t{b.end}\t"
                         f"{b.observed_total}\t{b.observed_top1}\n")
    print(f"wrote {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

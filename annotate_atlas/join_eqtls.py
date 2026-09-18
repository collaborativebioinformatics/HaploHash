#!/usr/bin/env python3
"""Count GTEx v8 significant eQTLs per haploblock, for AVI-ranking validation.

Hypothesis: the more constrained a block (high AVI top1_fraction), the fewer
eQTLs it should carry -- purifying selection removes regulatory variants
before they can reach the allele frequency GTEx needs to detect them.

Reads every "<Tissue>.v8.signif_variant_gene_pairs.txt.gz" file (GTEx v8
variant_id format: chr1_13550_G_A_b38, already hg38) and, per block, counts:
  eqtl_variant_count  -- unique variants significant in >=1 tissue
  eqtl_pair_count     -- unique (variant, gene) pairs, union across tissues
  eqtl_gene_count     -- unique eGenes
  eqtl_tissue_count   -- number of distinct tissues contributing >=1 hit

Usage:
  python3 join_eqtls.py blocks.bed GTEX_EQTL_DIR eqtl_by_block.tsv
"""

from __future__ import annotations

import bisect
import gzip
import glob
import os
import sys
from dataclasses import dataclass, field


@dataclass
class Block:
    chrom: str
    start: int
    end: int
    block_id: str
    variants: set = field(default_factory=set)
    pairs: set = field(default_factory=set)
    genes: set = field(default_factory=set)
    tissues: set = field(default_factory=set)


def load_blocks(bed_path: str) -> dict[str, list[Block]]:
    by_chrom: dict[str, list[Block]] = {}
    with open(bed_path) as fh:
        for line in fh:
            chrom, start, end, block_id = line.rstrip("\n").split("\t")
            by_chrom.setdefault(chrom, []).append(Block(chrom, int(start), int(end), block_id))
    for blocks in by_chrom.values():
        blocks.sort(key=lambda b: b.start)
    return by_chrom


def find_block(blocks: list[Block], starts: list[int], pos0: int) -> Block | None:
    i = bisect.bisect_right(starts, pos0) - 1
    if i < 0:
        return None
    block = blocks[i]
    return block if block.start <= pos0 < block.end else None


def main(argv: list[str]) -> int:
    if len(argv) != 4:
        print(f"usage: {argv[0]} blocks.bed GTEX_EQTL_DIR eqtl_by_block.tsv", file=sys.stderr)
        return 1
    bed_path, eqtl_dir, out_path = argv[1], argv[2], argv[3]

    by_chrom = load_blocks(bed_path)
    starts_by_chrom = {c: [b.start for b in blocks] for c, blocks in by_chrom.items()}

    pair_files = sorted(glob.glob(os.path.join(eqtl_dir, "*.v8.signif_variant_gene_pairs.txt.gz")))
    if not pair_files:
        print(f"no *.v8.signif_variant_gene_pairs.txt.gz under {eqtl_dir}", file=sys.stderr)
        return 1

    total_rows = 0
    unmatched = 0
    for path in pair_files:
        tissue = os.path.basename(path).split(".v8.signif_variant_gene_pairs.txt.gz")[0]
        with gzip.open(path, "rt") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            variant_idx = header.index("variant_id")
            gene_idx = header.index("gene_id")
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                variant_id, gene_id = fields[variant_idx], fields[gene_idx]
                chrom, pos_text = variant_id.split("_", 2)[:2]
                pos0 = int(pos_text) - 1
                total_rows += 1

                blocks = by_chrom.get(chrom)
                if blocks is None:
                    unmatched += 1
                    continue
                block = find_block(blocks, starts_by_chrom[chrom], pos0)
                if block is None:
                    unmatched += 1
                    continue

                block.variants.add(variant_id)
                block.pairs.add((variant_id, gene_id))
                block.genes.add(gene_id)
                block.tissues.add(tissue)
        print(f"[{tissue}] processed", file=sys.stderr)

    fieldnames = ["block_id", "chrom", "start", "end", "length_bp",
                  "eqtl_variant_count", "eqtl_pair_count", "eqtl_gene_count",
                  "eqtl_tissue_count", "eqtl_variants_per_kb"]
    with open(out_path, "w") as fh:
        fh.write("\t".join(fieldnames) + "\n")
        for blocks in by_chrom.values():
            for block in blocks:
                length_kb = (block.end - block.start) / 1000.0
                density = len(block.variants) / length_kb if length_kb > 0 else 0.0
                fh.write("\t".join(str(v) for v in [
                    block.block_id, block.chrom, block.start, block.end,
                    block.end - block.start, len(block.variants), len(block.pairs),
                    len(block.genes), len(block.tissues), round(density, 4),
                ]) + "\n")

    print(f"wrote {out_path}", file=sys.stderr)
    print(f"tissues processed: {len(pair_files)}", file=sys.stderr)
    print(f"total significant pairs read: {total_rows}, unmatched to any block: {unmatched}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

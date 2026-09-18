#!/usr/bin/env python3
"""Count GENCODE genes overlapping each haploblock -- a real gene-density
covariate, to test whether it (not selection) explains why high-AVI blocks
have more eQTLs.

Reads a GENCODE GTF (gzipped), keeps "gene" feature rows, and for every
haploblock in blocks.bed counts genes overlapping the block interval and
genes whose TSS (5' end, strand-aware) falls inside it.

Usage:
  python3 gene_density.py blocks.bed gencode.annotation.gtf.gz gene_density_by_block.tsv
"""

from __future__ import annotations

import bisect
import gzip
import sys


def load_blocks(bed_path: str) -> dict[str, list[list]]:
    by_chrom: dict[str, list[list]] = {}
    with open(bed_path) as fh:
        for line in fh:
            chrom, start, end, block_id = line.rstrip("\n").split("\t")
            by_chrom.setdefault(chrom, []).append([int(start), int(end), block_id, 0, 0])
    for rows in by_chrom.values():
        rows.sort()
    return by_chrom


def load_genes(gtf_path: str) -> dict[str, list[tuple[int, int, int]]]:
    """Returns chrom -> list of (start0, end, tss0) sorted by start, GTF 1-based inclusive -> BED half-open."""
    genes: dict[str, list[tuple[int, int, int]]] = {}
    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] != "gene":
                continue
            chrom = fields[0]
            start0, end = int(fields[3]) - 1, int(fields[4])
            strand = fields[6]
            tss0 = start0 if strand == "+" else end - 1
            genes.setdefault(chrom, []).append((start0, end, tss0))
    for rows in genes.values():
        rows.sort()
    return genes


def count_overlaps(blocks: dict[str, list[list]], genes: dict[str, list[tuple[int, int, int]]]) -> None:
    """Blocks are non-overlapping and sorted by start (checked upstream), so both
    the overlap window and the expiry window advance monotonically -- a single
    sweep per chromosome with a min-heap of active gene ends suffices."""
    import heapq

    for chrom, block_rows in blocks.items():
        gene_rows = genes.get(chrom, [])
        tss_values = sorted(g[2] for g in gene_rows)
        add_idx = 0
        active: list[int] = []  # heap of ends for genes with start < current block.end
        for block in block_rows:
            bstart, bend = block[0], block[1]
            while add_idx < len(gene_rows) and gene_rows[add_idx][0] < bend:
                heapq.heappush(active, gene_rows[add_idx][1])
                add_idx += 1
            while active and active[0] <= bstart:
                heapq.heappop(active)
            block[3] = len(active)
            lo_t = bisect.bisect_left(tss_values, bstart)
            hi_t = bisect.bisect_left(tss_values, bend)
            block[4] = hi_t - lo_t


def main(argv: list[str]) -> int:
    if len(argv) != 4:
        print(f"usage: {argv[0]} blocks.bed gencode.gtf.gz output.tsv", file=sys.stderr)
        return 1
    bed_path, gtf_path, out_path = argv[1:]

    blocks = load_blocks(bed_path)
    print(f"loaded {sum(len(v) for v in blocks.values())} blocks", file=sys.stderr)
    genes = load_genes(gtf_path)
    print(f"loaded {sum(len(v) for v in genes.values())} GENCODE genes", file=sys.stderr)
    count_overlaps(blocks, genes)

    with open(out_path, "w") as fh:
        fh.write("block_id\tchrom\tstart\tend\tn_genes_overlap\tn_genes_tss\n")
        for chrom, rows in blocks.items():
            for start, end, block_id, n_overlap, n_tss in sorted(rows, key=lambda r: r[2]):
                fh.write(f"{block_id}\t{chrom}\t{start}\t{end}\t{n_overlap}\t{n_tss}\n")
    total_genes = sum(r[3] for rows in blocks.values() for r in rows)
    print(f"wrote {out_path}: {total_genes} block-gene overlaps total", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

#!/usr/bin/env python3
"""Average recombination rate per haploblock, from the Eagle/Broad hg38
genetic map (Bherer et al. sex-averaged map).

For each block, linearly interpolates the genetic map's cumulative cM value
at block.start and block.end, then rate = (cM_end - cM_start) / length_Mb.
This is the standard way to get a regional recombination-rate estimate from
a sparse genetic map (average rate over the interval, not a point sample).

Usage:
  python3 recomb_rate.py blocks.bed genetic_map_hg38_withX.txt.gz recomb_rate_by_block.tsv
"""

from __future__ import annotations

import bisect
import gzip
import sys


def load_map(path: str) -> dict[str, tuple[list[int], list[float]]]:
    """Returns chrom -> (sorted positions, matching cumulative cM), chrom canonicalized to 'chrN'."""
    by_chrom: dict[str, list[tuple[int, float]]] = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        next(fh)  # header
        for line in fh:
            chrom, pos, _rate, cm = line.split()
            chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
            by_chrom.setdefault(chrom, []).append((int(pos), float(cm)))
    out = {}
    for chrom, rows in by_chrom.items():
        rows.sort()
        out[chrom] = ([p for p, _ in rows], [c for _, c in rows])
    return out


def interpolate(positions: list[int], cms: list[float], pos: int) -> float:
    if pos <= positions[0]:
        return cms[0]
    if pos >= positions[-1]:
        return cms[-1]
    i = bisect.bisect_left(positions, pos)
    if positions[i] == pos:
        return cms[i]
    x0, x1 = positions[i - 1], positions[i]
    y0, y1 = cms[i - 1], cms[i]
    return y0 + (y1 - y0) * (pos - x0) / (x1 - x0)


def main(argv: list[str]) -> int:
    if len(argv) != 4:
        print(f"usage: {argv[0]} blocks.bed genetic_map.txt.gz output.tsv", file=sys.stderr)
        return 1
    bed_path, map_path, out_path = argv[1:]

    gmap = load_map(map_path)
    print(f"loaded genetic map for {len(gmap)} chromosomes", file=sys.stderr)

    n_written = n_missing = 0
    with open(bed_path) as fh, open(out_path, "w") as out:
        out.write("block_id\tchrom\tstart\tend\trecomb_rate_cm_per_mb\n")
        for line in fh:
            chrom, start, end, block_id = line.rstrip("\n").split("\t")
            start, end = int(start), int(end)
            if chrom not in gmap:
                n_missing += 1
                continue
            positions, cms = gmap[chrom]
            cm_start = interpolate(positions, cms, start)
            cm_end = interpolate(positions, cms, end)
            length_mb = (end - start) / 1_000_000.0
            rate = (cm_end - cm_start) / length_mb if length_mb > 0 else 0.0
            out.write(f"{block_id}\t{chrom}\t{start}\t{end}\t{rate:.6g}\n")
            n_written += 1
    print(f"wrote {n_written} blocks ({n_missing} skipped, no map for chrom) to {out_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

#!/usr/bin/env python3
"""
haplohash.py - encode, decode and k-anonymity suppression for haploblock hashes.

Bit layout, 64-bit word, 53 used and 11 reserved:

    bits  field      range                      why
    ----  ---------  -------------------------  ---------------------------
     1    strand     hap0 / hap1                2 values
     5    chrom      1-24                       24 values
    16    block      0-65,535                   39,077 blocks observed
    13    cluster    0-8,191                    5,090 max cluster id observed
    10    variants   0-1,023                    as in the original scheme
     8    avi        quantised PHRED            256 levels
    11    reserved                              headroom

The original layout in the paper spends 4 bits on strand, 10 on chromosome and
20 each on block and cluster, which uses the whole word with nothing left over.
Tightening those to what the data actually needs frees 19 bits, which is where
the AVI score goes.

Encoding is reversible by construction: every field is a plain integer packed
at a known offset. That is the point of the attack half of this work. A hash
with a singleton cluster id identifies one haplotype, so suppression happens
before encoding, not after.

Usage:
    haplohash.py demo
    haplohash.py encode --clusters assignments.tsv --scores block_scores.tsv \\
                        --blocks block_stats.tsv --avi-col avi_top1_mean \\
                        --k 2 -o hashes.tsv
    haplohash.py decode --hashes hashes.tsv -o decoded.tsv
    haplohash.py suppress --blocks block_stats.tsv --k 5
"""

import argparse
import sys

# ----------------------------------------------------------------- bit layout
FIELDS = [
    ("strand",    1),
    ("chrom",     5),
    ("block",    16),
    ("cluster",  13),
    ("variants", 10),
    ("avi",       8),
]
RESERVED = 64 - sum(w for _, w in FIELDS)

_OFF, _o = {}, 0
for _name, _w in FIELDS:
    _OFF[_name] = (_o, _w)
    _o += _w


def pack(strand, chrom, block, cluster, variants, avi):
    """Pack six integers into one 64-bit word."""
    h = 0
    for name, value in (("strand", strand), ("chrom", chrom), ("block", block),
                        ("cluster", cluster), ("variants", variants), ("avi", avi)):
        off, width = _OFF[name]
        if not 0 <= value < (1 << width):
            raise ValueError(
                "{}={} does not fit in {} bits (max {})".format(
                    name, value, width, (1 << width) - 1))
        h |= value << off
    return h


def unpack(h):
    """Recover the six fields. This is the attack: no key, no secret."""
    out = {}
    for name, (off, width) in _OFF.items():
        out[name] = (h >> off) & ((1 << width) - 1)
    return out


# ----------------------------------------------------------------- quantising
def quantise(value, lo, hi, bits=8):
    """Linear quantisation of a float onto an integer range."""
    if hi <= lo:
        return 0
    n = (1 << bits) - 1
    q = round((value - lo) / (hi - lo) * n)
    return max(0, min(n, q))


def dequantise(q, lo, hi, bits=8):
    n = (1 << bits) - 1
    return lo + (q / n) * (hi - lo)


CHROM_MAP = {str(i): i for i in range(1, 23)}
CHROM_MAP.update({"X": 23, "Y": 24})


def chrom_to_int(c):
    c = str(c).replace("chr", "")
    if c not in CHROM_MAP:
        raise ValueError("unrecognised chromosome: " + str(c))
    return CHROM_MAP[c]


# ----------------------------------------------------------------- suppression
def suppression_report(block_stats_path, k):
    """
    How much data a k-anonymity threshold costs.

    k=2 (suppress singletons) is exact from singleton_count. Higher k needs the
    full cluster size distribution, which block_stats.tsv does not carry, so
    those are reported as lower bounds.
    """
    import pandas as pd
    df = pd.read_csv(block_stats_path, sep="\t")
    n_hap = int(df["n_haplotypes"].max())
    total_obs = n_hap * len(df)
    singles = int(df["singleton_count"].sum())

    print("blocks {:,}   haplotypes {:,}   observations {:,}".format(
        len(df), n_hap, total_obs))
    print()
    print("k = 2, suppress singleton clusters")
    print("  suppressed observations : {:,}  ({:.2f}% of data)".format(
        singles, 100 * singles / total_obs))
    print("  blocks affected         : {:,}  ({:.1f}%)".format(
        int((df["singleton_count"] > 0).sum()),
        100 * (df["singleton_count"] > 0).mean()))
    print()
    if k > 2:
        print("k = {}: exact cost needs the full cluster size distribution.".format(k))
        print("  lower bound is the k=2 figure above, since every cluster under")
        print("  size {} also includes every singleton.".format(k))
    print()
    print("for comparison, suppressing whole blocks until no singletons remain")
    n_blocks = int((df["singleton_count"] > 0).sum())
    print("  blocks removed : {:,} of {:,}  ({:.1f}%)".format(
        n_blocks, len(df), 100 * n_blocks / len(df)))


# ----------------------------------------------------------------- demo
def demo():
    print("bit layout, {} bits used, {} reserved\n".format(64 - RESERVED, RESERVED))
    for name, (off, width) in _OFF.items():
        print("  {:<9} offset {:>2}  width {:>2}  max {:,}".format(
            name, off, width, (1 << width) - 1))

    print("\nround trip")
    cases = [
        dict(strand=1, chrom=chrom_to_int("chr20"), block=38001,
             cluster=5090, variants=17, avi=quantise(31.4, 0, 60)),
        dict(strand=0, chrom=chrom_to_int("chr1"), block=0,
             cluster=0, variants=0, avi=0),
        dict(strand=1, chrom=chrom_to_int("chrX"), block=65535,
             cluster=8191, variants=1023, avi=255),
    ]
    ok = True
    for c in cases:
        h = pack(**c)
        back = unpack(h)
        match = all(back[k] == v for k, v in c.items())
        ok &= match
        print("  0x{:016x}  {}".format(h, "ok" if match else "MISMATCH"))
        if not match:
            print("    in ", c)
            print("    out", back)
    print("\nall round trips {}".format("passed" if ok else "FAILED"))

    print("\nthe attack, no key required")
    h = pack(strand=1, chrom=chrom_to_int("chr20"), block=38001,
             cluster=5090, variants=17, avi=quantise(31.4, 0, 60))
    f = unpack(h)
    print("  hash 0x{:016x}".format(h))
    print("  -> chrom {}  block {}  cluster {}  strand hap{}".format(
        f["chrom"], f["block"], f["cluster"], f["strand"]))
    print("  -> avi PHRED approx {:.1f}".format(dequantise(f["avi"], 0, 60)))
    print("  if cluster {} has one member, that is one named haplotype.".format(
        f["cluster"]))

    return 0 if ok else 1


# ----------------------------------------------------------------- cli
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd")

    sub.add_parser("demo", help="bit layout, round-trip test, decode example")

    s = sub.add_parser("suppress", help="k-anonymity cost report")
    s.add_argument("--blocks", required=True, help="block_stats.tsv")
    s.add_argument("--k", type=int, default=2)

    e = sub.add_parser("encode", help="encode cluster assignments to hashes")
    e.add_argument("--clusters", required=True,
                   help="TSV: block_id, haplotype_id, cluster_id, n_variants")
    e.add_argument("--scores", required=True, help="block_scores.tsv")
    e.add_argument("--blocks", required=True, help="block_stats.tsv")
    e.add_argument("--avi-col", default="avi_top1_mean",
                   help="which AVI column to pack (default avi_top1_mean)")
    e.add_argument("--k", type=int, default=2,
                   help="suppress clusters smaller than k (default 2)")
    e.add_argument("-o", "--out", required=True)

    d = sub.add_parser("decode", help="unpack hashes back to fields")
    d.add_argument("--hashes", required=True)
    d.add_argument("-o", "--out", required=True)

    a = p.parse_args()
    if a.cmd == "demo" or a.cmd is None:
        return demo()
    if a.cmd == "suppress":
        return suppression_report(a.blocks, a.k)
    if a.cmd == "encode":
        return encode_file(a)
    if a.cmd == "decode":
        return decode_file(a)


def encode_file(a):
    import pandas as pd
    cl = pd.read_csv(a.clusters, sep="\t")
    sc = pd.read_csv(a.scores, sep="\t")
    st = pd.read_csv(a.blocks, sep="\t")

    if a.avi_col not in sc.columns:
        sys.exit("column {} not in {}. available: {}".format(
            a.avi_col, a.scores, ", ".join(c for c in sc.columns if "avi" in c)))

    lo, hi = float(sc[a.avi_col].min()), float(sc[a.avi_col].max())
    print("quantising {} over [{:.3f}, {:.3f}] into 8 bits".format(a.avi_col, lo, hi))

    block_ids = {b: i for i, b in enumerate(sorted(st["block"].unique()))}
    avi = dict(zip(sc["block_id"], sc[a.avi_col]))

    sizes = cl.groupby(["block_id", "cluster_id"]).size()

    rows, dropped = [], 0
    for r in cl.itertuples(index=False):
        if sizes[(r.block_id, r.cluster_id)] < a.k:
            dropped += 1
            continue
        strand = 1 if str(r.haplotype_id).endswith("hap1") else 0
        chrom = chrom_to_int(str(r.block_id).split("_")[0])
        q = quantise(float(avi.get(r.block_id, lo)), lo, hi)
        h = pack(strand, chrom, block_ids[r.block_id], int(r.cluster_id),
                 int(getattr(r, "n_variants", 0)), q)
        rows.append((r.haplotype_id, r.block_id, "0x{:016x}".format(h)))

    with open(a.out, "w") as fh:
        fh.write("haplotype_id\tblock_id\thash\n")
        for x in rows:
            fh.write("\t".join(map(str, x)) + "\n")
    print("wrote {:,} hashes, suppressed {:,} at k={}".format(len(rows), dropped, a.k))


def decode_file(a):
    with open(a.hashes) as fh, open(a.out, "w") as out:
        header = fh.readline().rstrip("\n").split("\t")
        hi = header.index("hash")
        out.write("\t".join(header + list(_OFF)) + "\n")
        n = 0
        for line in fh:
            f = line.rstrip("\n").split("\t")
            d = unpack(int(f[hi], 16))
            out.write("\t".join(f + [str(d[k]) for k in _OFF]) + "\n")
            n += 1
    print("decoded {:,} hashes".format(n))


if __name__ == "__main__":
    sys.exit(main())

import csv
import tempfile
import unittest
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

import qc_variant_parquet as qc
from score_blocks import Block


def block_row(start=100, end=110):
    length = end - start
    return {
        "block_id": f"chr22_{start}_{end}", "chrom": "chr22",
        "start": start, "end": end, "length_bp": length,
        "n_positions": length, "n_variant_alleles": length * 3,
        "n_scored": length * 3, "n_missing_avi": 0, "avi_top1_count": 1,
    }


def variant(pos=101, start=100, end=110, **overrides):
    row = {
        "chrom": "chr22", "pos": pos, "ref": "A", "alt": "G",
        "variant_id": f"chr22-{pos}-A-G", "AC": 1, "AN": 5096, "AF": 0.0,
        "block_id": f"chr22_{start}-{end}", "block_start": start, "block_end": end,
        "avi_raw": 0.5, "avi_phred": 20.0,
    }
    row.update(overrides)
    return row


class PilotQCTest(unittest.TestCase):
    def setUp(self):
        self.blocks = qc.select_blocks(
            [block_row(), block_row(110, 120)],
            [Block("chr22", 100, 110, "first"), Block("chr22", 110, 120, "second")],
        )

    def test_filters_on_ac_and_preserves_reported_frequency(self):
        rows = [variant(), variant(102, AC=0), variant(103, AC=5095, AF=1.0)]
        observed, excluded, counts = qc.audit_variants(rows, self.blocks, snv=True)
        self.assertEqual([r["pos"] for r in observed], [101, 103])
        self.assertEqual(excluded[0]["exclusion_reason"], "AC_zero")
        self.assertEqual(observed[0]["AF"], 0.0)
        self.assertAlmostEqual(observed[0]["AF_from_counts"], 1 / 5096)
        self.assertAlmostEqual(observed[1]["MAF_from_counts"], 1 / 5096)
        self.assertEqual(counts[("chr22", 110, 120)]["observed_alleles"], 0)
        self.assertEqual(rows[0], variant())  # no mutation of the input

    def test_coordinate_join_ignores_block_id_separator(self):
        observed, _, _ = qc.audit_variants([variant(chrom="22")], self.blocks, snv=True)
        self.assertEqual(observed[0]["block_id"], "chr22_100-110")
        self.assertEqual(observed[0]["block_id_avi"], "chr22_100_110")
        self.assertEqual(observed[0]["block_key"], "chr22:100-110")

    def test_boundaries_use_one_based_positions_and_half_open_blocks(self):
        rows = [variant(110), variant(111, start=110, end=120)]
        observed, _, _ = qc.audit_variants(rows, self.blocks, snv=True)
        self.assertEqual([r["block_id_avi"] for r in observed], ["chr22_100_110", "chr22_110_120"])
        for pos in (100, 111):
            with self.subTest(pos=pos), self.assertRaisesRegex(ValueError, "outside"):
                qc.audit_variants([variant(pos)], self.blocks, snv=True)

    def test_deduplicates_alleles_without_collapsing_different_alts(self):
        rows = [variant(), variant(), variant(alt="T", variant_id="chr22-101-A-T")]
        observed, excluded, _ = qc.audit_variants(rows, self.blocks, snv=True)
        self.assertEqual(len(observed), 2)
        self.assertEqual({r["alt"] for r in observed}, {"G", "T"})
        self.assertEqual(excluded[0]["exclusion_reason"], "duplicate_allele")

    def test_conflicting_duplicate_fails(self):
        with self.assertRaisesRegex(ValueError, "conflicting duplicate"):
            qc.audit_variants([variant(), variant(AC=2)], self.blocks, snv=True)

    def test_missing_scores_are_retained_and_not_imputed_to_zero(self):
        rows = [variant(avi_phred=None), variant(102, avi_phred=float("nan"), avi_raw=None)]
        observed, _, counts = qc.audit_variants(rows, self.blocks, snv=True)
        self.assertEqual(counts[("chr22", 100, 110)]["observed_missing_avi"], 2)
        self.assertEqual(counts[("chr22", 100, 110)]["observed_scored"], 0)
        self.assertTrue(all(r["avi_phred"] is None and not r["has_avi_score"] for r in observed))

    def test_rejects_invalid_counts_alleles_and_scores(self):
        cases = [
            {"AC": -1}, {"AC": 5097}, {"AC": 1.5}, {"AN": 0},
            {"AF": float("nan")}, {"AF": 1.1}, {"alt": "G,T"},
            {"alt": "A"}, {"alt": "AG"}, {"avi_phred": float("inf")},
            {"avi_phred": -1}, {"pos": None},
        ]
        for case in cases:
            with self.subTest(case=case), self.assertRaises(ValueError):
                qc.audit_variants([variant(**case)], self.blocks, snv=True)

    def test_missing_and_duplicate_block_rows_fail(self):
        expected = [Block("chr22", 100, 110, "pilot")]
        for rows in ([], [block_row(), block_row()]):
            with self.subTest(rows=rows), self.assertRaises(ValueError):
                qc.select_blocks(rows, expected)
        inconsistent = {**block_row(), "n_scored": 31}
        with self.assertRaisesRegex(ValueError, "inconsistent AVI"):
            qc.select_blocks([inconsistent], expected)

    def test_indels_are_a_separate_audit(self):
        indel = {k: v for k, v in variant(alt="AG").items() if k in qc.INDEL_COLUMNS}
        observed, _, _ = qc.audit_variants([indel], self.blocks, snv=False)
        self.assertEqual(len(observed), 1)
        self.assertNotIn("has_avi_score", observed[0])
        with self.assertRaises(ValueError):
            qc.audit_variants([indel], self.blocks, snv=True)

    def test_cli_outputs_and_no_outputs_on_failed_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            blocks = qc.pilot_blocks()
            first = blocks[0]
            snvs = root / "snvs.parquet"
            indels = root / "indels.parquet"
            scores = root / "scores.tsv"
            rows = [
                variant(first.start + 1, first.start, first.end),
                variant(first.start + 2, first.start, first.end, AC=0),
            ]
            pq.write_table(pa.Table.from_pylist(rows), snvs)
            indel = {k: v for k, v in rows[0].items() if k in qc.INDEL_COLUMNS}
            indel["alt"] = "AG"
            pq.write_table(pa.Table.from_pylist([indel]), indels)
            block_rows = [block_row(b.start, b.end) for b in blocks]
            qc.write_tsv(scores, block_rows, list(block_rows[0]))
            out = root / "qc"
            args = ["--snvs", str(snvs), "--indels", str(indels), "--block-scores", str(scores), "--outdir", str(out)]
            self.assertEqual(qc.main(args), 0)
            cleaned = pq.read_table(out / "chr22_observed_snvs.parquet").to_pylist()
            self.assertEqual(len(cleaned), 1)
            self.assertEqual(cleaned[0]["AF"], 0.0)
            self.assertAlmostEqual(cleaned[0]["AF_from_counts"], 1 / 5096)
            self.assertEqual(
                {p.name for p in out.iterdir()},
                {"chr22_observed_snvs.parquet", "chr22_excluded_snvs.tsv", "chr22_block_qc.tsv"},
            )
            with (out / "chr22_block_qc.tsv").open() as stream:
                summaries = list(csv.DictReader(stream, delimiter="\t"))
            self.assertEqual(len(summaries), 12)
            self.assertEqual(sum(int(r["snv_observed_alleles"]) for r in summaries), 1)
            self.assertEqual(sum(int(r["snv_AC_zero"]) for r in summaries), 1)
            rows[0]["pos"] = first.start
            pq.write_table(pa.Table.from_pylist(rows), snvs)
            bad_out = root / "failed"
            with self.assertRaisesRegex(ValueError, "outside"):
                qc.main([*args[:-1], str(bad_out)])
            self.assertFalse(bad_out.exists())


if __name__ == "__main__":
    unittest.main()

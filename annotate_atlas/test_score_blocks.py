import argparse
import tempfile
import unittest
from pathlib import Path

import score_blocks


class ScoreBlocksTest(unittest.TestCase):
    def test_resolves_named_header_columns(self):
        args = argparse.Namespace(
            chrom_column=None, position_column=None, phred_column=None
        )
        self.assertEqual(
            score_blocks.resolve_columns(
                args, ["chromosome", "position", "alternate", "AVI_PHRED"]
            ),
            (0, 1, 3),
        )

    def test_headerless_input_requires_explicit_phred_column(self):
        args = argparse.Namespace(
            chrom_column=None, position_column=None, phred_column=None
        )
        with self.assertRaisesRegex(ValueError, "could not infer the phred column"):
            score_blocks.resolve_columns(args, None)

    def test_aggregates_alleles_with_half_open_bed_boundaries(self):
        blocks = [
            score_blocks.Block("chr1", 0, 10, "b1"),
            score_blocks.Block("chr1", 10, 20, "b2"),
        ]
        rows = [
            "chr1\t1\tA\t21",
            "chr1\t1\tC\t10",
            "chr1\t10\tG\t30",
            "chr1\t11\tT\t25",
            "chr1\t12\tA\t.",
        ]

        summaries, unmatched = score_blocks.aggregate_rows(
            rows,
            blocks,
            chrom_column=0,
            position_column=1,
            phred_column=3,
            top_k=2,
        )

        self.assertEqual(unmatched, 0)
        first, second = summaries
        self.assertEqual(first.n_positions, 2)
        self.assertEqual(first.n_variant_alleles, 3)
        self.assertEqual(first.n_top1, 2)
        self.assertEqual(first.maximum, 30)
        self.assertEqual(sorted(first._top_scores), [21, 30])
        self.assertAlmostEqual(first.sum_all / first.n_scored, (21 + 10 + 30) / 3)
        self.assertEqual(second.n_positions, 2)
        self.assertEqual(second.n_variant_alleles, 2)
        self.assertEqual(second.n_scored, 1)
        self.assertEqual(second.n_missing_avi, 1)

    def test_load_blocks_rejects_overlap(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "blocks.bed"
            path.write_text("chr1\t0\t10\tb1\nchr1\t9\t20\tb2\n")
            with self.assertRaisesRegex(ValueError, "overlapping blocks"):
                score_blocks.load_blocks(path)


if __name__ == "__main__":
    unittest.main()

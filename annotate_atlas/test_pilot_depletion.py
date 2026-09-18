import tempfile
import unittest
from pathlib import Path

import pilot_depletion as depletion
from qc_variant_parquet import select_blocks
from score_blocks import Block
from test_qc_variant_parquet import block_row, variant


class DepletionTest(unittest.TestCase):
    def blocks(self, high_counts=(3, 6)):
        return select_blocks(
            [{**block_row(start, start + 10), "avi_top1_count": high}
             for start, high in zip((100, 110), high_counts)],
            [Block("chr22", 100, 110, "a"), Block("chr22", 110, 120, "b")],
        )

    def test_block_expectations_and_pooled_ratio(self):
        variants = [variant(), variant(102, avi_phred=19.99),
                    variant(111, 110, 120, avi_phred=30)]
        rows, frequencies = depletion.summarize(variants, self.blocks())
        # Two scored alleles at 10% opportunity; one at 20% opportunity.
        self.assertEqual([r["observed_top1"] for r in rows], [1, 1])
        self.assertEqual([r["expected_top1"] for r in rows], [0.2, 0.2])
        self.assertEqual(depletion.pooled_totals(rows), (2, 0.4, 5.0))
        self.assertEqual(len(frequencies["high"]), 2)
        self.assertEqual(len(frequencies["low"]), 1)
        # The pooled ratio must use summed expectations, not mean block ratios.
        totals = depletion.pooled_totals([
            {"observed_top1": 1, "expected_top1": 2},
            {"observed_top1": 3, "expected_top1": 1},
        ])
        self.assertAlmostEqual(totals[2], 4 / 3)

    def test_missing_scores_do_not_enter_either_side_of_comparison(self):
        rows, frequencies = depletion.summarize(
            [variant(), variant(102, avi_phred=None)], self.blocks())
        self.assertEqual(rows[0]["n_observed"], 2)
        self.assertEqual(rows[0]["n_missing_avi"], 1)
        self.assertEqual(rows[0]["n_observed_scored"], 1)
        self.assertEqual(rows[0]["expected_top1"], 0.1)
        self.assertEqual(frequencies["low"], [])

    def test_zero_opportunity_and_empty_blocks_have_undefined_ratios(self):
        rows, _ = depletion.summarize([variant(avi_phred=0)], self.blocks((0, 6)))
        self.assertEqual([r["expected_top1"] for r in rows], [0.0, 0.0])
        self.assertTrue(all(r["observed_expected_ratio"] is None for r in rows))
        self.assertIsNone(depletion.pooled_totals(rows)[2])

    def test_unscored_block_has_no_expectation(self):
        blocks = self.blocks()
        blocks[("chr22", 100, 110)].update(n_scored=0, n_missing_avi=30, avi_top1_count=0)
        rows, _ = depletion.summarize([variant(avi_phred=None)], blocks)
        self.assertIsNone(rows[0]["expected_top1"])

    def test_frequencies_use_counts_and_handle_fixed_alternate_alleles(self):
        rows, frequencies = depletion.summarize(
            [variant(AC=5095, AF=1.0), variant(102, AC=5096, AF=1.0)], self.blocks())
        self.assertEqual(rows[0]["observed_top1"], 2)
        self.assertAlmostEqual(frequencies["high"][0], 1 / 5096)
        self.assertEqual(frequencies["high"][1], 0)
        with tempfile.TemporaryDirectory() as directory:
            depletion.plot_results(rows, frequencies, Path(directory))
            self.assertEqual(len(list(Path(directory).glob("*.png"))), 2)

    def test_rejects_unclean_inputs_and_impossible_counts(self):
        for variants, blocks in [
            ([variant(AC=0)], self.blocks()),
            ([variant(), variant()], self.blocks()),
            ([variant()], self.blocks((0, 6))),
        ]:
            with self.subTest(variants=variants), self.assertRaises(ValueError):
                depletion.summarize(variants, blocks)


if __name__ == "__main__":
    unittest.main()

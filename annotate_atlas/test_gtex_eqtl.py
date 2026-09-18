import gzip
import io
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

import gtex_eqtl as g


class GTExTest(unittest.TestCase):
    def test_official_calls_override_rounded_cutoffs_and_require_tested_coverage(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            payload = gzip.compress(b"variant_id\nchr1_101_A_G_b38\n")
            with tarfile.open(root / "GTEx_Analysis_v8_eQTL.tar", "w") as archive:
                info = tarfile.TarInfo("GTEx_Analysis_v8_eQTL/test.v8.signif_variant_gene_pairs.txt.gz")
                info.size = len(payload)
                archive.addfile(info, io.BytesIO(payload))
            variants = pd.DataFrame({"variant": ["chr1_101_A_G", "chr2_201_C_T"],
                                     "is_eqtl": [False, True], "is_lead": [False, False]})
            with patch.object(g, "DATA", root):
                g.validate_hits("test", variants)
                self.assertEqual(variants.is_eqtl.tolist(), [True, False])
                with self.assertRaisesRegex(ValueError, "absent from tested background"):
                    g.validate_hits("test", variants.iloc[1:].copy())

    def test_pair_collapse_deduplicates_genes_and_accepts_missing_allele_counts(self):
        genes = pd.DataFrame({"gene_id": ["A", "B"], "tss": [100, 150],
                              "egene": [True, False], "pval_nominal_threshold": [0.01, 0.01],
                              "lead_variant": ["chr1_101_A_G", "chr1_101_A_G"]})
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "allpairs").mkdir()
            with gzip.open(root / "allpairs/test.tsv.gz", "wt") as stream:
                stream.write("variant\tchromosome\tposition\tmaf\tgene_id\tpvalue\tan\ttype\tref\talt\n")
                stream.write("chr1_101_A_G\t1\t101\t0.1\tA\t0.001\t100\tSNP\tA\tG\n" * 2)
                stream.write("chr1_101_A_G\t1\t101\t0.1\tB\t0.2\tNaN\tSNP\tA\tG\n")
                stream.write("chr2_201_C_T\t2\t201\t0.2\tA\t0.1\t100\tSNP\tC\tT\n")
                stream.write("chr1_102_A_AT\t1\t102\t0.1\tA\t0.001\t100\tINDEL\tA\tAT\n")
            target = root / "variants.parquet"
            with patch.multiple(g, WORK=root, DATA=root), patch.object(g, "load_genes", return_value=genes):
                g.collapse_variants("test", target)
            variants = pd.read_parquet(target).set_index("variant")
            self.assertEqual(len(variants), 2)
            row = variants.loc["chr1_101_A_G"]
            self.assertEqual(row.n_tested_genes, 2)
            self.assertEqual(row.min_tss_distance, 1)
            self.assertEqual(row.an, 100)
            self.assertTrue(row.is_eqtl and row.is_lead)
            self.assertFalse(variants.loc["chr2_201_C_T", "is_eqtl"])

    def test_interval_assignment_keeps_boundary_and_gap_conventions(self):
        blocks = pd.DataFrame({"chrom": ["chr1", "chr1", "chr2"], "start": [100, 110, 0],
                               "end": [110, 120, 10], "block_index": [0, 1, 2]})
        variants = pd.DataFrame({"chrom": ["chr1"] * 4 + ["chr2", "chr3"],
                                 "position": [100, 101, 110, 111, 10, 10]})
        np.testing.assert_array_equal(g.assign_blocks(variants, blocks), [-1, 0, 0, 1, 2, -1])

    def test_tissue_aliases_must_resolve_uniquely(self):
        self.assertEqual(g.canonical_tissue("Skin_Not_Sun_Exposed", ["Skin_Not_Sun_Exposed_Suprapubic"]),
                         "Skin_Not_Sun_Exposed_Suprapubic")
        with self.assertRaises(ValueError):
            g.canonical_tissue("Brain", ["Brain_Cortex", "Brain_Cerebellum"])

    def test_adjustment_removes_a_pure_covariate_composition_effect(self):
        blocks = pd.DataFrame({"block_index": np.arange(10), "avi_quintile": np.tile(np.arange(1, 6), 2),
                               "region": [f"chr1:{i}" for i in range(10)]})
        rows = []
        for block in blocks.itertuples():
            high_proportion = block.avi_quintile / 6
            for stratum, tested, rate in [(0, round(600 * (1 - high_proportion)), 0.1),
                                           (1, round(600 * high_proportion), 0.9)]:
                rows.append(dict(block_index=block.block_index, stratum=stratum, tested=tested,
                                 eqtl=round(tested * rate), lead=round(tested * rate), finemapped=round(tested * rate)))
        _, quintiles, summary = g.analyze_strata(pd.DataFrame(rows), blocks, "test", n_bootstrap=100)
        np.testing.assert_allclose(quintiles.observed_expected, 1)
        np.testing.assert_allclose(summary.adjusted_high_low_ratio, 1)
        self.assertTrue((summary.raw_high_low_ratio > 3).all())

    def test_zero_case_outcomes_do_not_create_spurious_enrichment(self):
        blocks = pd.DataFrame({"block_index": np.arange(5), "avi_quintile": np.arange(1, 6),
                               "region": [f"chr1:{i}" for i in range(5)]})
        strata = pd.DataFrame({"block_index": np.arange(5), "stratum": 0, "tested": 100,
                               "eqtl": 0, "lead": 0, "finemapped": 0})
        _, _, summary = g.analyze_strata(strata, blocks, "test", n_bootstrap=100)
        self.assertTrue(summary.adjusted_high_low_ratio.isna().all())
        self.assertTrue(summary.ci_low.isna().all())


if __name__ == "__main__":
    unittest.main()

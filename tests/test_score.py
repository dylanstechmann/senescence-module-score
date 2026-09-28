import io
import json
import os
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from senescore.cli import bakeoff, main as cli_main
from senescore.genes import FRIDMAN_SENESCENCE, GENE_SETS, ORTHOGONAL, SASP_COPPE, SENMAYO_HUMAN, get_gene_set
from senescore.io import read_expression_csv
from senescore.score import (
    ScoreError,
    module_score,
    pearson_correlation,
    spearman_correlation,
)
from senescore.synthetic import spiked_cohort


class ScoreTests(unittest.TestCase):
    def test_gene_set_size_and_orthogonal_markers_are_not_members(self):
        self.assertEqual(len(SENMAYO_HUMAN), 125)
        self.assertEqual(len(set(SENMAYO_HUMAN)), 125)
        for gene in ORTHOGONAL:
            self.assertNotIn(gene, SENMAYO_HUMAN)

    def test_alternative_gene_sets_in_registry(self):
        self.assertIn("senmayo", GENE_SETS)
        self.assertIn("fridman", GENE_SETS)
        self.assertIn("sasp", GENE_SETS)
        self.assertEqual(len(FRIDMAN_SENESCENCE), len(set(FRIDMAN_SENESCENCE)))
        self.assertEqual(len(SASP_COPPE), len(set(SASP_COPPE)))
        self.assertEqual(get_gene_set("fridman")["symbols"], FRIDMAN_SENESCENCE)
        self.assertEqual(get_gene_set("sasp")["symbols"], SASP_COPPE)
        with self.assertRaises(KeyError):
            get_gene_set("unknown_set")

    def test_spike_beats_a_random_gene_set(self):
        report = bakeoff(0)
        self.assertGreater(report["senmayo_cohen_d"], 2.0)
        self.assertLess(abs(report["random_set_cohen_d"]), 1.0)
        self.assertGreater(report["senmayo_cohen_d"], abs(report["random_set_cohen_d"]) + 1.5)

    def test_correlation_metrics(self):
        a = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        b = np.array([2.0, 4.0, 6.0, 8.0, 10.0])
        self.assertAlmostEqual(pearson_correlation(a, b), 1.0)
        self.assertAlmostEqual(spearman_correlation(a, b), 1.0)
        c = np.array([5.0, 4.0, 3.0, 2.0, 1.0])
        self.assertAlmostEqual(pearson_correlation(a, c), -1.0)
        self.assertAlmostEqual(spearman_correlation(a, c), -1.0)

    def test_low_coverage_raises(self):
        matrix, genes, _ = spiked_cohort(seed=1)
        keep = [gene for gene in genes if gene not in SENMAYO_HUMAN[:80]]
        cols = [genes.index(gene) for gene in keep]
        with self.assertRaises(ScoreError):
            module_score(matrix[:, cols], keep, SENMAYO_HUMAN)

    def test_csv_roundtrip_ranks_the_spiked_sample_last(self):
        matrix, genes, labels = spiked_cohort(n_per_group=6, seed=2)
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "expr.csv")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write("sample," + ",".join(genes) + "\n")
                for i, row in enumerate(matrix):
                    vals = ",".join(f"{v:.4f}" for v in row)
                    handle.write(f"s{i},{vals}\n")
            ids, read_genes, read_matrix = read_expression_csv(path)
        self.assertEqual(ids[0], "s0")
        result = module_score(read_matrix, read_genes, SENMAYO_HUMAN, seed=2)
        spiked_scores = result["scores"][labels == "spiked"]
        reference_scores = result["scores"][labels == "reference"]
        self.assertGreater(float(spiked_scores.mean()), float(reference_scores.mean()))
        self.assertIn("CDKN1A", result["orthogonal_z"])
        # orthogonal markers were not spiked, so their group gap stays small
        z = np.array(result["orthogonal_z"]["CDKN1A"])
        self.assertLess(abs(float(z[labels == "spiked"].mean() - z[labels == "reference"].mean())), 1.0)

    def test_cli_compare_and_gene_set_scoring(self):
        rng = np.random.default_rng(3)
        all_genes = sorted(list(set(SENMAYO_HUMAN) | set(FRIDMAN_SENESCENCE) | {f"BG{i}" for i in range(200)}))
        n_samples = 10
        matrix = rng.uniform(1.0, 5.0, size=(n_samples, len(all_genes)))
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "expr.csv")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write("sample," + ",".join(all_genes) + "\n")
                for i, row in enumerate(matrix):
                    vals = ",".join(f"{v:.4f}" for v in row)
                    handle.write(f"sample_{i},{vals}\n")

            # test score with --gene-set fridman
            stdout_score = io.StringIO()
            with patch("sys.stdout", stdout_score):
                rc = cli_main(["score", path, "--gene-set", "fridman"])
                self.assertEqual(rc, 0)
            score_data = json.loads(stdout_score.getvalue())
            self.assertEqual(score_data["gene_set"], "fridman")
            self.assertIn("fridman_module_score", score_data["samples"][0])

            # test compare subcommand
            stdout_cmp = io.StringIO()
            with patch("sys.stdout", stdout_cmp):
                rc = cli_main(["compare", path, "--set-a", "senmayo", "--set-b", "fridman"])
                self.assertEqual(rc, 0)
            cmp_data = json.loads(stdout_cmp.getvalue())
            self.assertIn("comparison", cmp_data)
            self.assertIn("correlation", cmp_data)
            self.assertEqual(cmp_data["comparison"]["set_a"]["key"], "senmayo")
            self.assertEqual(cmp_data["comparison"]["set_b"]["key"], "fridman")
            self.assertIn("pearson_r", cmp_data["correlation"])
            self.assertIn("spearman_rho", cmp_data["correlation"])
            self.assertEqual(len(cmp_data["samples"]), len(matrix))


if __name__ == "__main__":
    unittest.main()


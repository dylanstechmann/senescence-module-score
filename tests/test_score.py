import io
import hashlib
import json
import os
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from senescore.cli import bakeoff, main as cli_main
from senescore.genes import (
    FRIDMAN_SENESCENCE,
    FRIDMAN_SENESCENCE_DOWN,
    FRIDMAN_SENESCENCE_UP,
    GENE_SETS,
    ORTHOGONAL,
    SASP_COPPE,
    SENMAYO_HUMAN,
    get_gene_set,
)
from senescore.io import read_expression_csv
from senescore.score import (
    ScoreError,
    module_score,
    pearson_correlation,
    spearman_correlation,
    signed_direction_score,
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
        self.assertEqual(get_gene_set("fridman")["source_status"], "custom_unverified")
        self.assertEqual(get_gene_set("sasp")["source_status"], "custom_unverified")
        self.assertEqual(get_gene_set("senmayo")["source_status"], "published_gene_set")
        self.assertEqual(len(FRIDMAN_SENESCENCE_UP), 77)
        self.assertEqual(len(FRIDMAN_SENESCENCE_DOWN), 13)
        self.assertEqual(get_gene_set("fridman_up")["systematic_id"], "M9143")
        self.assertEqual(get_gene_set("fridman_down")["systematic_id"], "M9487")
        self.assertIn("COL3A1", FRIDMAN_SENESCENCE_DOWN)
        self.assertIn("EGR1", FRIDMAN_SENESCENCE_DOWN)
        self.assertNotIn("COL3A1", FRIDMAN_SENESCENCE_UP)
        self.assertNotIn("EGR1", FRIDMAN_SENESCENCE_UP)
        self.assertEqual(get_gene_set("fridman_up")["source_status"], "published_gene_set")
        self.assertEqual(get_gene_set("fridman_down")["source_status"], "published_gene_set")
        with self.assertRaises(KeyError):
            get_gene_set("unknown_set")

    def test_signed_direction_score_subtracts_down_from_up(self):
        np.testing.assert_array_equal(signed_direction_score([3, -1], [1, 2]), [2, -3])
        for up, down in [([1], []), ([1, np.nan], [1, 2]), ([[1]], [[1]])]:
            with self.subTest(up=up, down=down), self.assertRaises(ScoreError):
                signed_direction_score(up, down)

    def test_synthetic_quiescence_and_proliferation_spikes_do_not_mimic_signature_spike(self):
        rng = np.random.default_rng(881)
        quiescence = [f"QUIESCENCE_FIXTURE_{i}" for i in range(8)]
        proliferation = [f"PROLIFERATION_FIXTURE_{i}" for i in range(8)]
        genes = list(SENMAYO_HUMAN) + quiescence + proliferation + [f"BG{i}" for i in range(300)]
        gene_means = rng.uniform(2.0, 4.0, size=len(genes))
        matrix = gene_means[None, :] + rng.normal(0.0, 0.08, size=(40, len(genes)))
        index = {gene: i for i, gene in enumerate(genes)}
        groups = {
            "reference": np.arange(0, 10),
            "quiescence_only": np.arange(10, 20),
            "proliferation_only": np.arange(20, 30),
            "senmayo_spike": np.arange(30, 40),
        }
        for rows, fixture_signature in ((groups["quiescence_only"], quiescence),
                                        (groups["proliferation_only"], proliferation),
                                        (groups["senmayo_spike"], SENMAYO_HUMAN)):
            matrix[np.ix_(rows, [index[g] for g in fixture_signature])] += 2.0
        background_positions = [index[f"BG{i}"] for i in range(300)]
        matrix[np.ix_(groups["quiescence_only"], background_positions)] += 0.1
        matrix[np.ix_(groups["proliferation_only"], background_positions)] -= 0.1
        result = module_score(
            matrix,
            genes,
            SENMAYO_HUMAN,
            seed=5,
            n_ctrl=5,
            n_bins=10,
            block_from_controls=set(SENMAYO_HUMAN) | set(quiescence) | set(proliferation),
        )
        means = {name: float(result["scores"][rows].mean()) for name, rows in groups.items()}
        self.assertGreater(means["senmayo_spike"], means["reference"] + 1.0)
        self.assertLess(abs(means["quiescence_only"] - means["reference"]), 0.5)
        self.assertLess(abs(means["proliferation_only"] - means["reference"]), 0.5)

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

    def test_correlations_use_average_tie_ranks_and_report_undefined(self):
        self.assertAlmostEqual(spearman_correlation(np.array([1, 1, 2, 2]), np.array([1, 2, 1, 2])), 0.0)
        self.assertIsNone(spearman_correlation(np.ones(4), np.full(4, 3.0)))
        self.assertIsNone(pearson_correlation(np.ones(4), np.full(4, 3.0)))
        self.assertIsNone(spearman_correlation(np.array([1.0, np.nan]), np.array([1.0, 2.0])))

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
            self.assertEqual(score_data["gene_set_source_status"], "custom_unverified")
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
            self.assertEqual(cmp_data["comparison"]["set_a"]["source_status"], "published_gene_set")
            self.assertEqual(cmp_data["comparison"]["set_b"]["source_status"], "custom_unverified")
            self.assertIn("pearson_r", cmp_data["correlation"])
            self.assertIn("spearman_rho", cmp_data["correlation"])
            self.assertEqual(len(cmp_data["samples"]), len(matrix))

    def test_cli_hash_matches_bytes_used_when_input_changes_during_scoring(self):
        from senescore import cli

        matrix, genes, _ = spiked_cohort(n_per_group=2, seed=19)
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "expr.csv")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write("sample," + ",".join(genes) + "\n")
                for i, row in enumerate(matrix):
                    handle.write(f"s{i}," + ",".join(map(str, row)) + "\n")
            with open(path, "rb") as handle:
                parsed_bytes = handle.read()
            replacement = "\n".join(["sample,G", "a,999", "b,1000", "c,1001", "d,1002"]) + "\n"
            original_score = cli.module_score

            def replace_after_parse(*args, **kwargs):
                result = original_score(*args, **kwargs)
                with open(path, "w", encoding="utf-8") as handle:
                    handle.write(replacement)
                return result

            stdout = io.StringIO()
            with patch("senescore.cli.module_score", side_effect=replace_after_parse), patch("sys.stdout", stdout):
                self.assertEqual(cli_main(["score", path]), 0)
            payload = json.loads(stdout.getvalue())
            self.assertEqual(payload["input_sha256"], hashlib.sha256(parsed_bytes).hexdigest())
            self.assertNotEqual(payload["input_sha256"], hashlib.sha256(replacement.encode()).hexdigest())


if __name__ == "__main__":
    unittest.main()

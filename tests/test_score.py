import os
import tempfile
import unittest

import numpy as np

from senescore.cli import bakeoff
from senescore.genes import ORTHOGONAL, SENMAYO_HUMAN
from senescore.io import read_expression_csv
from senescore.score import ScoreError, module_score
from senescore.synthetic import spiked_cohort


class ScoreTests(unittest.TestCase):
    def test_gene_set_size_and_orthogonal_markers_are_not_members(self):
        self.assertEqual(len(SENMAYO_HUMAN), 125)
        self.assertEqual(len(set(SENMAYO_HUMAN)), 125)
        for gene in ORTHOGONAL:
            self.assertNotIn(gene, SENMAYO_HUMAN)

    def test_spike_beats_a_random_gene_set(self):
        report = bakeoff(0)
        self.assertGreater(report["senmayo_cohen_d"], 2.0)
        self.assertLess(abs(report["random_set_cohen_d"]), 1.0)
        self.assertGreater(report["senmayo_cohen_d"], abs(report["random_set_cohen_d"]) + 1.5)

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


if __name__ == "__main__":
    unittest.main()

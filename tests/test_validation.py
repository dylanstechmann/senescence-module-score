import tempfile
import unittest
from pathlib import Path

import numpy as np

from senescore.io import read_expression_csv
from senescore.score import ScoreError, cohen_d, module_score


class ValidationTests(unittest.TestCase):
    def test_rejects_nonfinite_duplicates_and_empty_signature(self):
        cases = [(np.array([[np.nan, 1]]), ["A", "B"], ["A"]),
                 (np.ones((2, 2)), ["A", " a "], ["A"]),
                 (np.ones((2, 2)), ["A", "B"], [])]
        for matrix, genes, signature in cases:
            with self.assertRaises(ScoreError):
                module_score(matrix, genes, signature)

    def test_control_provenance_reconstructs_score(self):
        matrix = np.array([[1, 3, 2, 4], [2, 4, 1, 3]], dtype=float)
        genes = ["A", "B", "C", "D"]
        result = module_score(matrix, genes, [" a "], n_bins=1, n_ctrl=2)
        controls = result["control_genes"]["A"]
        self.assertNotIn("A", controls)
        expected = matrix[:, 0] - matrix[:, [genes.index(g) for g in controls]].mean(axis=1)
        np.testing.assert_allclose(result["scores"], expected)
        order = [3, 1, 0, 2]
        reordered = module_score(matrix[:, order], [genes[i] for i in order], ["A"], n_bins=1, n_ctrl=2)
        np.testing.assert_allclose(result["scores"], reordered["scores"])

    def test_csv_errors_are_actionable(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.csv"
            for text in ["sample,A,B\ns1,1\ns2,2,3\n", "sample,A,B\ns1,1,2\ns1,3,4\n",
                         "sample,A,B\ns1,nan,2\ns2,3,4\n"]:
                path.write_text(text)
                with self.assertRaises(ScoreError):
                    read_expression_csv(path)

    def test_zero_variance_is_not_reported_as_no_effect(self):
        with self.assertRaisesRegex(ScoreError, "zero pooled variance"):
            cohen_d(np.array([1, 1]), np.array([2, 2]))

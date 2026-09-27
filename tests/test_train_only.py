import unittest

import numpy as np

from senescore.genes import SENMAYO_HUMAN
from senescore.score import ScoreError, module_score, module_score_train_only
from senescore.synthetic import spiked_cohort


class TrainOnlyTests(unittest.TestCase):
    def test_all_rows_match_the_cohort_score(self):
        matrix, genes, _labels = spiked_cohort(seed=0)
        full = module_score(matrix, genes, SENMAYO_HUMAN, seed=0, block_from_controls=set(SENMAYO_HUMAN))
        held = module_score_train_only(
            matrix, genes, SENMAYO_HUMAN, train_rows=np.arange(len(matrix)),
            seed=0, block_from_controls=set(SENMAYO_HUMAN),
        )
        self.assertTrue(np.allclose(full["scores"], held["scores"]))
        self.assertEqual(full["control_genes"], held["control_genes"])

    def test_test_only_expression_does_not_change_training_controls(self):
        matrix, genes, _labels = spiked_cohort(n_per_group=8, seed=3)
        train = np.arange(0, 8)
        # Move one background gene only in the held-out rows so cohort-wide bins can change.
        background = next(i for i, gene in enumerate(genes) if gene not in SENMAYO_HUMAN)
        leaked = matrix.copy()
        leaked[8:, background] += 50.0
        cohort = module_score(leaked, genes, SENMAYO_HUMAN, seed=1, n_ctrl=5, n_bins=10,
                              block_from_controls=set(SENMAYO_HUMAN))
        frozen = module_score_train_only(
            leaked, genes, SENMAYO_HUMAN, train_rows=train, seed=1, n_ctrl=5, n_bins=10,
            block_from_controls=set(SENMAYO_HUMAN),
        )
        untouched = module_score_train_only(
            matrix, genes, SENMAYO_HUMAN, train_rows=train, seed=1, n_ctrl=5, n_bins=10,
            block_from_controls=set(SENMAYO_HUMAN),
        )
        self.assertEqual(frozen["control_genes"], untouched["control_genes"])
        self.assertNotEqual(cohort["control_genes"], frozen["control_genes"])

    def test_train_rows_are_validated(self):
        matrix, genes, _labels = spiked_cohort(seed=0)
        with self.assertRaises(ScoreError):
            module_score_train_only(matrix, genes, SENMAYO_HUMAN, train_rows=[0])
        with self.assertRaises(ScoreError):
            module_score_train_only(matrix, genes, SENMAYO_HUMAN, train_rows=[0, 0])


if __name__ == "__main__":
    unittest.main()

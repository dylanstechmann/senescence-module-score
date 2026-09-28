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
        self.assertEqual(full["orthogonal_z"], held["orthogonal_z"])

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

    def test_held_out_marker_values_do_not_set_training_z_scores(self):
        matrix, genes, _labels = spiked_cohort(n_per_group=8, seed=3)
        train = np.arange(0, 8)
        marker = genes.index("CDKN1A")
        changed = matrix.copy()
        changed[8:, marker] += 50.0
        frozen = module_score_train_only(matrix, genes, SENMAYO_HUMAN, train_rows=train)
        held_out_changed = module_score_train_only(changed, genes, SENMAYO_HUMAN, train_rows=train)
        np.testing.assert_allclose(
            frozen["orthogonal_z"]["CDKN1A"][:8],
            held_out_changed["orthogonal_z"]["CDKN1A"][:8],
        )
        self.assertEqual(held_out_changed["orthogonal_fit_on"], "train_rows_only")
        leaked = module_score(changed, genes, SENMAYO_HUMAN)
        self.assertFalse(np.allclose(leaked["orthogonal_z"]["CDKN1A"][:8],
                                         held_out_changed["orthogonal_z"]["CDKN1A"][:8]))

    def test_train_rows_are_validated(self):
        matrix, genes, _labels = spiked_cohort(seed=0)
        with self.assertRaises(ScoreError):
            module_score_train_only(matrix, genes, SENMAYO_HUMAN, train_rows=[0])
        with self.assertRaises(ScoreError):
            module_score_train_only(matrix, genes, SENMAYO_HUMAN, train_rows=[0, 0])
        for invalid in ([0.5, 1.5], [True, False]):
            with self.assertRaises(ScoreError):
                module_score_train_only(matrix, genes, SENMAYO_HUMAN, train_rows=invalid)

    def test_train_only_validation_and_gene_reordering(self):
        matrix, genes, _labels = spiked_cohort(n_per_group=8, seed=5)
        train = np.arange(8)
        first = module_score_train_only(matrix, genes, SENMAYO_HUMAN, train_rows=train)
        order = np.random.default_rng(7).permutation(len(genes))
        reordered = module_score_train_only(
            matrix[:, order], [genes[i] for i in order], SENMAYO_HUMAN, train_rows=train
        )
        np.testing.assert_allclose(first["scores"], reordered["scores"])
        self.assertEqual(first["control_genes"], reordered["control_genes"])
        with self.assertRaises(ScoreError):
            module_score_train_only(matrix, genes[:-1] + [genes[0]], SENMAYO_HUMAN, train_rows=train)
        invalid = matrix.copy()
        invalid[0, 0] = np.nan
        with self.assertRaises(ScoreError):
            module_score_train_only(invalid, genes, SENMAYO_HUMAN, train_rows=train)


if __name__ == "__main__":
    unittest.main()

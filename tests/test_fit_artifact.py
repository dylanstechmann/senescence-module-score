import csv
import hashlib
import io
import json
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

import numpy as np

from senescore.cli import main
from senescore.genes import SENMAYO_HUMAN
from senescore.score import (
    ScoreError,
    module_score_train_only,
    transform_module_score,
    validate_fit_artifact,
)
from senescore.synthetic import spiked_cohort


def _write_expression(path, sample_ids, genes, matrix):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample", *genes])
        for sample_id, values in zip(sample_ids, matrix):
            writer.writerow([sample_id, *[repr(float(value)) for value in values]])


class FitArtifactTests(unittest.TestCase):
    def setUp(self):
        self.matrix, self.genes, _labels = spiked_cohort(n_per_group=8, seed=3)
        self.sample_ids = [f"sample-{i}" for i in range(len(self.matrix))]
        self.train_rows = np.arange(8)

    def test_saved_fit_transforms_later_rows_without_refitting(self):
        fitted = module_score_train_only(
            self.matrix,
            self.genes,
            SENMAYO_HUMAN,
            train_rows=self.train_rows,
            seed=1,
            n_ctrl=5,
            n_bins=10,
            training_sample_ids=self.sample_ids[:8],
            gene_set_key="senmayo",
        )
        artifact = fitted["fit_artifact"]
        transformed = transform_module_score(self.matrix[8:], self.genes, artifact)
        np.testing.assert_allclose(transformed["scores"], fitted["scores"][8:])
        for marker in fitted["orthogonal_z"]:
            np.testing.assert_allclose(
                transformed["orthogonal_z"][marker], fitted["orthogonal_z"][marker][8:]
            )
        self.assertEqual(artifact["training_sample_ids"], self.sample_ids[:8])
        self.assertEqual(artifact["control_genes"], fitted["control_genes"])
        self.assertEqual(transformed["controls_fit_on"], "fit_artifact_training_rows")
        self.assertEqual(artifact["feature_schema_sha256"], hashlib.sha256(
            json.dumps(self.genes, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest())

    def test_held_out_changes_do_not_change_the_fitted_artifact(self):
        first = module_score_train_only(
            self.matrix, self.genes, SENMAYO_HUMAN, train_rows=self.train_rows, seed=4, n_bins=10
        )["fit_artifact"]
        changed = self.matrix.copy()
        changed[8:, 0] += 500
        second = module_score_train_only(
            changed, self.genes, SENMAYO_HUMAN, train_rows=self.train_rows, seed=4, n_bins=10
        )["fit_artifact"]
        self.assertEqual(first, second)

    def test_transform_fails_closed_on_incompatible_feature_schema(self):
        artifact = module_score_train_only(
            self.matrix, self.genes, SENMAYO_HUMAN, train_rows=self.train_rows
        )["fit_artifact"]
        order = np.arange(len(self.genes))[::-1]
        with self.assertRaisesRegex(ScoreError, "feature schema/order"):
            transform_module_score(self.matrix[:, order], [self.genes[i] for i in order], artifact)
        broken = dict(artifact)
        broken["control_genes"] = dict(artifact["control_genes"])
        gene = next(iter(broken["control_genes"]))
        broken["control_genes"][gene] = ["NOT_A_FEATURE"]
        with self.assertRaisesRegex(ScoreError, "invalid controls"):
            validate_fit_artifact(broken)

    def test_cli_exports_and_reuses_fit_artifact(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            training_csv = root / "training.csv"
            later_csv = root / "later.csv"
            training_ids = root / "training-ids.txt"
            artifact_path = root / "fit.json"
            _write_expression(training_csv, self.sample_ids, self.genes, self.matrix)
            _write_expression(later_csv, self.sample_ids[8:], self.genes, self.matrix[8:])
            training_ids.write_text("\n".join(self.sample_ids[:8]) + "\n", encoding="utf-8")

            fitted_output = io.StringIO()
            with redirect_stdout(fitted_output):
                self.assertEqual(main([
                    "score", str(training_csv), "--train-samples", str(training_ids),
                    "--fit-artifact", str(artifact_path), "--bins", "10",
                ]), 0)
            fit_payload = json.loads(fitted_output.getvalue())
            self.assertTrue(artifact_path.exists())
            artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
            self.assertEqual(artifact["training_sample_ids"], self.sample_ids[:8])
            self.assertEqual(fit_payload["training_sample_ids_sha256"], artifact["training_sample_ids_sha256"])

            applied_output = io.StringIO()
            with redirect_stdout(applied_output):
                self.assertEqual(main([
                    "score", str(later_csv), "--apply-artifact", str(artifact_path),
                ]), 0)
            applied = json.loads(applied_output.getvalue())
            expected = [row["senmayo_module_score"] for row in fit_payload["samples"][8:]]
            actual = [row["senmayo_module_score"] for row in applied["samples"]]
            self.assertEqual(actual, expected)
            self.assertEqual(applied["fit_artifact_version"], 1)


if __name__ == "__main__":
    unittest.main()

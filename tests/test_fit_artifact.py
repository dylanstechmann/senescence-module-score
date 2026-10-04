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
from senescore.genes import CITATION, SENMAYO_HUMAN
from senescore.score import (
    ScoreError,
    fit_module_score,
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

    def test_registered_gene_set_provenance_is_bound_to_signature_and_citation(self):
        artifact = fit_module_score(
            self.matrix, self.genes, SENMAYO_HUMAN, self.train_rows,
            n_bins=10, gene_set_key="SeNMayo", gene_set_source_status="published_gene_set",
            citation=CITATION,
        )
        self.assertEqual(artifact["gene_set_key"], "senmayo")
        self.assertEqual(artifact["gene_set_source_status"], "published_gene_set")
        validate_fit_artifact(artifact)
        with self.assertRaisesRegex(ScoreError, "does not match its registered gene set"):
            fit_module_score(self.matrix, self.genes, SENMAYO_HUMAN[:-1], self.train_rows,
                             n_bins=10, gene_set_key="senmayo", citation=CITATION)
        with self.assertRaisesRegex(ScoreError, "citation conflicts"):
            fit_module_score(self.matrix, self.genes, SENMAYO_HUMAN, self.train_rows,
                             n_bins=10, gene_set_key="senmayo", citation="unrelated citation")

    def test_published_status_requires_registry_key_and_registered_status_cannot_be_custom(self):
        with self.assertRaisesRegex(ScoreError, "requires a registered gene_set_key"):
            fit_module_score(self.matrix, self.genes, SENMAYO_HUMAN, self.train_rows,
                             n_bins=10, gene_set_source_status="published_gene_set", citation=CITATION)
        with self.assertRaisesRegex(ScoreError, "conflicts with its registered"):
            fit_module_score(self.matrix, self.genes, SENMAYO_HUMAN, self.train_rows,
                             n_bins=10, gene_set_key="senmayo", gene_set_source_status="custom_unverified",
                             citation=CITATION)

    def test_unregistered_panel_is_labeled_custom_unverified_and_artifact_rejects_tampering(self):
        signature = self.genes[:12]
        artifact = fit_module_score(
            self.matrix, self.genes, signature, self.train_rows, n_bins=10,
            gene_set_key="local-panel-v1", citation="Local panel; source not verified.",
        )
        self.assertEqual(artifact["gene_set_source_status"], "custom_unverified")
        validate_fit_artifact(artifact)
        for field, value, message in (
            ("gene_set_source_status", "published_gene_set", "unregistered gene sets"),
            ("gene_set_key", None, "published gene-set status requires"),
        ):
            broken = dict(artifact)
            broken[field] = value
            if field == "gene_set_key":
                broken["gene_set_source_status"] = "published_gene_set"
            with self.subTest(field=field), self.assertRaisesRegex(ScoreError, message):
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
            self.assertEqual(artifact["gene_set_source_status"], "published_gene_set")
            self.assertEqual(fit_payload["gene_set_source_status"], "published_gene_set")
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
            self.assertEqual(applied["gene_set"], "senmayo")
            self.assertEqual(applied["gene_set_source_status"], "published_gene_set")


if __name__ == "__main__":
    unittest.main()

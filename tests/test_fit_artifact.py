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
from senescore.genes import CITATION, FRIDMAN_SENESCENCE_DOWN, FRIDMAN_SENESCENCE_UP, SENMAYO_HUMAN
from senescore.score import (
    ScoreError,
    fit_module_score,
    fit_signed_fridman_score,
    module_score_train_only,
    transform_module_score,
    transform_signed_fridman_score,
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

    def test_signed_fridman_artifact_keeps_directional_fit_and_is_train_only(self):
        rng = np.random.default_rng(19)
        genes = list(FRIDMAN_SENESCENCE_UP + FRIDMAN_SENESCENCE_DOWN) + [f"BG{i}" for i in range(180)]
        matrix = rng.uniform(1.0, 5.0, size=(12, len(genes)))
        # Test data carry a clear UP and DOWN directional shift; no claims about biology.
        index = {gene: i for i, gene in enumerate(genes)}
        matrix[8:, [index[gene] for gene in FRIDMAN_SENESCENCE_UP]] += 2.0
        matrix[8:, [index[gene] for gene in FRIDMAN_SENESCENCE_DOWN]] -= 1.0
        fit = fit_signed_fridman_score(
            matrix, genes, np.arange(8), seed=5, n_ctrl=3, n_bins=10,
            training_sample_ids=self.sample_ids[:8], training_input_sha256="a" * 64,
        )
        result = transform_signed_fridman_score(matrix, genes, fit)
        np.testing.assert_allclose(result["scores"], result["up_scores"] - result["down_scores"])
        self.assertEqual(result["up_coverage"], 1.0)
        self.assertEqual(result["down_coverage"], 1.0)
        self.assertEqual(fit["up_fit"]["gene_set_key"], "fridman_up")
        self.assertEqual(fit["down_fit"]["gene_set_key"], "fridman_down")

        changed = matrix.copy()
        changed[8:, :] *= 100
        second = fit_signed_fridman_score(
            changed, genes, np.arange(8), seed=5, n_ctrl=3, n_bins=10,
            training_sample_ids=self.sample_ids[:8], training_input_sha256="a" * 64,
        )
        self.assertEqual(fit, second)

        broken = json.loads(json.dumps(fit))
        broken["up_fit"]["signature"][0] = "COL3A1"
        with self.assertRaisesRegex(ScoreError, "invalid fridman_up component"):
            transform_signed_fridman_score(matrix, genes, broken)

    def test_signed_fit_requires_60_percent_coverage_on_each_direction(self):
        rng = np.random.default_rng(324)
        for up, down in ((FRIDMAN_SENESCENCE_UP[:45], FRIDMAN_SENESCENCE_DOWN),
                         (FRIDMAN_SENESCENCE_UP, FRIDMAN_SENESCENCE_DOWN[:7])):
            genes = list(up + down) + [f"BG{i}" for i in range(200)]
            matrix = rng.uniform(1, 4, size=(10, len(genes)))
            with self.subTest(up=len(up), down=len(down)), self.assertRaisesRegex(ScoreError, "signature is in the table"):
                fit_signed_fridman_score(matrix, genes, list(range(6)), n_bins=10)

    def test_signed_fit_rejects_mixed_training_references_and_signature_controls(self):
        rng = np.random.default_rng(75)
        genes = list(FRIDMAN_SENESCENCE_UP + FRIDMAN_SENESCENCE_DOWN) + [f"BG{i}" for i in range(150)]
        matrix = rng.uniform(1, 4, size=(10, len(genes)))
        fitted = fit_signed_fridman_score(matrix, genes, list(range(6)), n_bins=10)
        # Metadata is optional for Python callers: a fit without a CSV hash can apply.
        transform_signed_fridman_score(matrix, genes, fitted)
        other = fit_signed_fridman_score(matrix, genes, list(range(2, 8)), n_bins=10)
        mixed = json.loads(json.dumps(fitted))
        mixed["down_fit"] = other["down_fit"]
        with self.assertRaisesRegex(ScoreError, "components disagree"):
            transform_signed_fridman_score(matrix, genes, mixed)
        broken = json.loads(json.dumps(fitted))
        broken["training_input_sha256"] = "b" * 64
        with self.assertRaisesRegex(ScoreError, "wrapper disagrees"):
            transform_signed_fridman_score(matrix, genes, broken)
        broken = json.loads(json.dumps(fitted))
        first = next(iter(broken["up_fit"]["control_genes"]))
        broken["up_fit"]["control_genes"][first] = [FRIDMAN_SENESCENCE_DOWN[0]]
        with self.assertRaisesRegex(ScoreError, "directional signature genes"):
            transform_signed_fridman_score(matrix, genes, broken)

    def test_fit_artifact_rejects_false_coverage_missing_metadata_and_invalid_controls(self):
        artifact = fit_module_score(self.matrix, self.genes, SENMAYO_HUMAN, self.train_rows, n_bins=10)
        for field, value, expected in (
            ("coverage", 0.2, "coverage"),
            ("missing", ["FAKE"], "missing-gene"),
        ):
            broken = json.loads(json.dumps(artifact))
            broken[field] = value
            with self.subTest(field=field), self.assertRaisesRegex(ScoreError, expected):
                transform_module_score(self.matrix, self.genes, broken)
        gene = next(iter(artifact["control_genes"]))
        for controls in ([gene], [artifact["control_genes"][gene][0]] * 2):
            broken = json.loads(json.dumps(artifact))
            broken["control_genes"][gene] = controls
            with self.subTest(controls=controls), self.assertRaisesRegex(ScoreError, "invalid controls"):
                validate_fit_artifact(broken)

    def test_compare_training_scores_and_controls_ignore_held_out_changes(self):
        rng = np.random.default_rng(633)
        signature = list(dict.fromkeys(SENMAYO_HUMAN + FRIDMAN_SENESCENCE_UP))
        genes = signature + [f"BG{i}" for i in range(350)]
        matrix = rng.uniform(1, 4, size=(10, len(genes)))
        ids = [f"s-{i}" for i in range(10)]
        with tempfile.TemporaryDirectory() as temp:
            expression = Path(temp) / "expression.csv"
            sample_file = Path(temp) / "train.txt"
            sample_file.write_text("\n".join(ids[:6]), encoding="utf-8")
            def compare(values):
                _write_expression(expression, ids, genes, values)
                output = io.StringIO()
                with redirect_stdout(output):
                    main(["compare", str(expression), "--train-samples", str(sample_file), "--bins", "10"])
                return json.loads(output.getvalue())
            first = compare(matrix)
            changed = matrix.copy()
            changed[6:] *= 100
            second = compare(changed)
        self.assertEqual(first["comparison"]["set_b"]["key"], "fridman_up")
        self.assertEqual(first["controls_fit_on"], "train_rows_only")
        self.assertEqual(first["training_sample_ids"], ids[:6])
        self.assertEqual(first["control_genes"], second["control_genes"])
        self.assertEqual(first["samples"][:6], second["samples"][:6])

    def test_signed_fridman_cli_round_trip(self):
        rng = np.random.default_rng(27)
        genes = list(FRIDMAN_SENESCENCE_UP + FRIDMAN_SENESCENCE_DOWN) + [f"BG{i}" for i in range(150)]
        matrix = rng.uniform(1.0, 4.0, size=(8, len(genes)))
        with tempfile.TemporaryDirectory() as temp:
            expression = Path(temp) / "expression.csv"
            sample_file = Path(temp) / "train.txt"
            artifact_file = Path(temp) / "fit.json"
            _write_expression(expression, [f"sample-{i}" for i in range(8)], genes, matrix)
            sample_file.write_text("\n".join(f"sample-{i}" for i in range(6)), encoding="utf-8")
            first = io.StringIO()
            with redirect_stdout(first):
                self.assertEqual(main(["fridman-signed", str(expression), "--train-samples", str(sample_file),
                                       "--fit-artifact", str(artifact_file), "--controls", "3", "--bins", "10"]), 0)
            initial = json.loads(first.getvalue())
            second = io.StringIO()
            with redirect_stdout(second):
                self.assertEqual(main(["fridman-signed", str(expression), "--apply-artifact", str(artifact_file)]), 0)
            applied = json.loads(second.getvalue())
            self.assertEqual(initial["gene_set_source_status"], "published_gene_set")
            self.assertEqual(initial["samples"], applied["samples"])
            self.assertEqual(initial["training_sample_ids"], [f"sample-{i}" for i in range(6)])
            self.assertEqual(initial["feature_schema_sha256"], applied["feature_schema_sha256"])
            self.assertEqual(applied["fit_artifact_hash_basis"], "file_bytes")
            self.assertEqual(applied["fit_artifact_sha256"], hashlib.sha256(artifact_file.read_bytes()).hexdigest())

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

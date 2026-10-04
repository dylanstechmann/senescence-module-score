"""Meaningful failure/leakage checks for the real GEO importer and evaluation."""

import gzip
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from senescore.genes import GENE_SETS, ORTHOGONAL
from senescore.public_validation import (
    descriptive_contrast, fit_panels, parse_geo_samples, read_counts, reduce_hgnc, sha256,
)
from senescore.score import ScoreError
from senescore.score import transform_module_score, transform_signed_fridman_score


ROOT = Path(__file__).resolve().parents[1]
SNAPSHOT = ROOT / "validation" / "GSE268487"
RAW = ROOT / "data" / "public" / "GSE268487" / "counts.txt.gz"


def metadata(accession="GSM1", condition="proliferating", replicate=1):
    treatment = {"proliferating": "None", "quiescent": "Contact Inhibition + Serum Starvation",
                 "senescent": "Replicative Exhaustion"}[condition]
    return (
        f"^SAMPLE = {accession}\n!Sample_title = LF1 cells, {condition}, replicate {replicate}\n"
        "!Sample_organism_ch1 = Homo sapiens\n!Sample_characteristics_ch1 = cell line: LF1\n"
        "!Sample_characteristics_ch1 = cell type: lung fibroblasts\n"
        f"!Sample_characteristics_ch1 = treatment: {treatment}\n"
    )


class PublicValidationTests(unittest.TestCase):
    def setUp(self):
        self.samples = parse_geo_samples(metadata() + metadata("GSM2", "quiescent"))
        self.mapping = {"ENSG00000000001": "A", "ENSG00000000003": "A", "ENSG00000000004": "ZERO"}

    def import_table(self, body):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "counts.tsv"
            path.write_text("\tqui1\tpro1\n" + body, encoding="utf-8")
            return read_counts(path, self.mapping, self.samples)

    def test_unmapped_gene_denominator_and_symbol_sum_before_log(self):
        matrix, raw, genes, samples, statistics = self.import_table(
            "ENSG00000000001.1\t20\t10\nENSG00000000002.3\t60\t70\n"
            "ENSG00000000003.1\t20\t20\nENSG00000000004.1\t0\t0\n"
            "__no_feature\t1000\t1000\ngSpikein_ERCC-001\t100\t100\n123\t50\t50\n"
        )
        self.assertEqual([s["accession"] for s in samples], ["GSM2", "GSM1"])
        self.assertEqual(genes, ["A", "ZERO"])
        np.testing.assert_array_equal(raw, [[40, 0], [30, 0]])
        np.testing.assert_allclose(matrix[:, 0], np.log2([400001, 300001]))
        self.assertEqual(statistics["ensembl_library_counts"], [100, 100])
        self.assertEqual(statistics["excluded_non_ensembl_counts"], [1150, 1150])
        self.assertEqual(statistics["unmapped_ensembl_rows"], 1)

    def test_fractional_negative_nan_and_infinite_counts_refused_including_excluded_rows(self):
        for value in ("-1", "1.2", "NaN", "inf", "broken"):
            with self.subTest(value=value), self.assertRaises(ScoreError):
                self.import_table(f"ENSG00000000001.1\t10\t10\n__no_feature\t{value}\t2\n")

    def test_duplicate_identifiers_and_version_variants_refused(self):
        for identifier in ("ENSG00000000001.1", "ENSG00000000001.2"):
            with self.subTest(identifier=identifier), self.assertRaisesRegex(ScoreError, "duplicate"):
                self.import_table(f"ENSG00000000001.1\t10\t10\n{identifier}\t5\t5\n")

    def test_pseudoautosomal_y_rows_are_valid_ensembl_loci_and_merge_by_symbol(self):
        matrix, raw, genes, samples, statistics = self.import_table(
            "ENSG00000000001.1\t20\t10\nENSG00000000001.1_PAR_Y\t30\t40\n"
            "ENSG00000000002.1\t50\t50\n"
        )
        self.assertEqual(statistics["ensembl_rows"], 3)
        self.assertEqual(statistics["ensembl_library_counts"], [100, 100])
        np.testing.assert_array_equal(raw, [[50], [50]])
        np.testing.assert_allclose(matrix[:, 0], np.log2([500001, 500001]))
        with self.assertRaisesRegex(ScoreError, "duplicate"):
            self.import_table("ENSG00000000001.1_PAR_Y\t10\t10\nENSG00000000001.2_PAR_Y\t20\t20\n")

    def test_ragged_rows_zero_libraries_and_malformed_identifiers_refused(self):
        for body in ("ENSG00000000001.1\t1\n", "ENSG00000000001.1\t0\t10\n", "ENSGbogus\t10\t10\n"):
            with self.subTest(body=body), self.assertRaises(ScoreError):
                self.import_table(body)

    def test_columns_must_match_actual_metadata_exactly(self):
        self.samples[0]["column"] = "invented"
        with self.assertRaisesRegex(ScoreError, "columns"):
            self.import_table("ENSG00000000001.1\t1\t2\n")

    def test_metadata_uses_titles_without_inventing_donor_or_pair_ids(self):
        samples = parse_geo_samples(metadata("GSM1", "senescent", 3))
        self.assertEqual(samples[0]["column"], "sen3")
        self.assertIsNone(samples[0]["donor_id"])
        self.assertIsNone(samples[0]["experiment_block"])
        self.assertEqual(samples[0]["reported_replicate"], 3)

    def test_metadata_retains_exact_geo_extraction_and_lane_processing_fields(self):
        protocol = (
            "Total RNA for proliferating and RS samples (extracted using trizol) and cell pellets for quiescent cells "
            "were sent to Azenta. where libraries were prepared using a stranded ribosomal depletion method. "
            "The proliferating and RS cells were split across two HiSeq lanes. The quiescent cells were run on one HiSeq lane."
        )
        library = "NEBNext® Ultra™ II Directional RNA Library Prep Kit for Illumina®"
        text = metadata() + f"!Sample_extract_protocol_ch1 = {protocol}\n!Sample_extract_protocol_ch1 = {library}\n"
        sample = parse_geo_samples(text)[0]
        self.assertEqual(sample["source_fields"]["!Sample_extract_protocol_ch1"], [protocol, library])
        committed = json.loads((SNAPSHOT / "samples.json").read_text())
        for record in committed:
            self.assertEqual(record["source_fields"]["!Sample_extract_protocol_ch1"], [protocol, library])

    def test_metadata_refuses_wrong_treatment_titles_and_duplicate_replicates(self):
        for text in (metadata().replace("treatment: None", "treatment: Wrong"),
                     metadata().replace("replicate 1", "replicate 4"),
                     metadata() + metadata("GSM2")):
            with self.subTest(text=text), self.assertRaises(ScoreError):
                parse_geo_samples(text)

    def test_hgnc_mapping_excludes_withdrawn_and_ambiguous_ids(self):
        mapping, ambiguous = reduce_hgnc(
            "ensembl_gene_id\tsymbol\tstatus\nENSG00000000001\tA\tApproved\n"
            "ENSG00000000001\tB\tApproved\nENSG00000000002\tC\tEntry Withdrawn\n"
            "ENSG00000000003\tD\tApproved\n"
        )
        self.assertEqual(mapping, {"ENSG00000000003": "D"})
        self.assertEqual(ambiguous, ["ENSG00000000001"])
        with self.assertRaises(ScoreError):
            reduce_hgnc("symbol\tstatus\nA\tApproved\n")

    def test_exact_small_sample_permutation_and_tied_auroc(self):
        result = descriptive_contrast([10, 11, 12, 1, 2, 3], [0, 1, 2], [3, 4, 5], permutation=True)
        self.assertEqual(result["directional_auroc"], 1)
        self.assertEqual(result["permutation_assignments"], 20)
        self.assertEqual(result["exact_two_sided_permutation_p"], .1)
        tied = descriptive_contrast([1, 1, 1, 1], [0, 1], [2, 3], permutation=True)
        self.assertEqual(tied["directional_auroc"], .5)
        self.assertEqual(tied["exact_two_sided_permutation_p"], 1)

    def test_held_out_expression_never_selects_panel_controls_or_marker_reference(self):
        names = sorted(set(ORTHOGONAL).union(*(set(entry["symbols"]) for entry in GENE_SETS.values())))
        names += [f"CONTROL_{i}" for i in range(800)]
        rng = np.random.default_rng(33)
        matrix = rng.normal(size=(9, len(names)))
        samples = [{"accession": f"GSM{i}", "condition": condition}
                   for i, condition in enumerate(["quiescent"] * 3 + ["proliferating"] * 3 + ["senescent"] * 3)]
        fits, refusals, train = fit_panels(matrix, names, samples, "a" * 64)
        changed = matrix.copy()
        changed[[0, 1, 2, 6, 7, 8]] += rng.normal(1000, 2000, (6, len(names)))
        later_fits, later_refusals, later_train = fit_panels(changed, names, samples, "a" * 64)
        self.assertEqual(fits, later_fits)
        self.assertEqual(refusals, later_refusals)
        self.assertEqual(train, later_train)
        self.assertEqual(fits["senmayo"]["training_sample_ids"], ["GSM3", "GSM4", "GSM5"])
        self.assertEqual(fits["fridman_signed"]["up_fit"]["training_sample_ids"], ["GSM3", "GSM4", "GSM5"])

    def test_committed_source_snapshots_match_pinned_hashes_and_have_real_accessions(self):
        manifest = json.loads((SNAPSHOT / "sources.json").read_text())
        for name, expected in manifest["snapshots"].items():
            self.assertEqual(sha256(SNAPSHOT / name), expected)
        samples = json.loads((SNAPSHOT / "samples.json").read_text())
        self.assertEqual(len(samples), 9)
        self.assertEqual([s["accession"] for s in samples], [f"GSM82918{i}" for i in range(59, 68)])
        self.assertTrue(all(s["donor_id"] is None and s["experiment_block"] is None for s in samples))
        recorded = json.loads((SNAPSHOT / "results.json").read_text())
        self.assertEqual(recorded["software"]["importer_source_sha256"], sha256(ROOT / "src" / "senescore" / "public_validation.py"))
        self.assertEqual(recorded["software"]["importer_correction_receipt_sha256"], sha256(SNAPSHOT / "IMPORTER_CORRECTION.md"))

    @unittest.skipUnless(RAW.exists(), "actual GEO matrix is intentionally outside git; run documented download for integration check")
    def test_actual_geo_import_with_pinned_hash_and_measured_signature_coverage(self):
        manifest = json.loads((SNAPSHOT / "sources.json").read_text())
        self.assertEqual(sha256(RAW), manifest["counts"]["sha256"])
        samples = json.loads((SNAPSHOT / "samples.json").read_text())
        full_metadata = RAW.parent / "metadata.soft"
        if full_metadata.exists():
            self.assertEqual(sha256(full_metadata), manifest["geo_metadata"]["sha256"])
            self.assertEqual(parse_geo_samples(full_metadata.read_text()), samples)
        with gzip.open(SNAPSHOT / "hgnc_mapping.tsv.gz", "rt") as handle:
            import csv
            mapping = {row["ensembl_gene_id"]: row["symbol"] for row in csv.DictReader(handle, delimiter="\t")}
        matrix, raw, names, ordered, statistics = read_counts(RAW, mapping, samples)
        self.assertEqual(statistics["ensembl_rows"], 58780)
        self.assertEqual(matrix.shape[0], 9)
        self.assertTrue(np.isfinite(matrix).all())
        self.assertEqual([s["column"] for s in ordered], ["qui1", "qui2", "qui3", "pro1", "pro2", "pro3", "sen1", "sen2", "sen3"])
        self.assertGreaterEqual(len(set(names) & set(GENE_SETS["senmayo"]["symbols"])), 75)
        fits, refusals, train = fit_panels(matrix, names, ordered, manifest["counts"]["sha256"])
        recorded = json.loads((SNAPSHOT / "results.json").read_text())
        self.assertEqual(refusals, recorded["refusals"])
        self.assertEqual([ordered[i]["accession"] for i in train], recorded["training_sample_ids"])
        for key, fit in fits.items():
            transform = transform_signed_fridman_score if key == "fridman_signed" else transform_module_score
            actual = transform(matrix, names, fit)["scores"]
            expected = [recorded["panels"][key]["scores_by_accession"][sample["accession"]] for sample in ordered]
            np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)


if __name__ == "__main__":
    unittest.main()

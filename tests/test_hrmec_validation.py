"""Checks for the prespecified GSE160356 evaluator, run on a SYNTHETIC fixture.

Nothing here touches the real GEO data. The fixture is a made-up nine-library table with the real
file layout (Ensembl identifiers, HTSeq summary rows, per-library gzip files in one tar) so the
importer, the leave-one-out baseline protocol, the verdict rule and the receipt can be checked
without calculating a single score on the real series.
"""

import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import shutil
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from senescore import hrmec_validation as hv
from senescore.genes import GENE_SETS, ORTHOGONAL, SENMAYO_HUMAN
from senescore.score import ScoreError, fit_module_score

ROOT = Path(__file__).resolve().parents[1]
GSMS = {"EP_CT": ["GSM9000001", "GSM9000004", "GSM9000007"],
        "LP_CT": ["GSM9000002", "GSM9000005", "GSM9000008"],
        "ETO": ["GSM9000003", "GSM9000006", "GSM9000009"]}
TREATMENT = {"EP_CT": "Early Passage", "LP_CT": "Late Passage", "ETO": "Etoposide Treatment"}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def sample_block(prefix, replicate, accession, *, title=None, treatment=None, organism="Homo sapiens",
                 member=None):
    member = member or f"{accession}_run{replicate}.genename.htcounts.txt.gz"
    return (
        f"^SAMPLE = {accession}\n!Sample_title = {title or f'{prefix}_{replicate}'}\n"
        "!Sample_source_name_ch1 = Human Retinal Microvascular Endothelial Cells\n"
        f"!Sample_organism_ch1 = {organism}\n"
        f"!Sample_characteristics_ch1 = treatment: {treatment or TREATMENT[prefix]}\n"
        "!Sample_characteristics_ch1 = tissue: Endothelial Cells\n"
        "!Sample_characteristics_ch1 = passage: P3 to P8\n"
        f"!Sample_supplementary_file_1 = ftp://example.invalid/{accession}/suppl/{member}\n")


def soft_text(blocks=None):
    series = ("^SERIES = GSE160356\n!Series_title = Synthetic fixture\n!Series_overall_design = Synthetic\n"
              "!Series_status = Public\n")
    if blocks is None:
        blocks = [sample_block(prefix, i + 1, accession) for prefix, items in GSMS.items()
                  for i, accession in enumerate(items)]
    return series + "".join(blocks)


def build_fixture(directory, *, shift=1.0, n_fillers=2200, seed=0):
    """Write a synthetic raw tar, snapshot directory and HGNC mapping; return their paths."""
    registered = sorted(set(SENMAYO_HUMAN).union(*(set(e["symbols"]) for e in GENE_SETS.values()), ORTHOGONAL))
    symbols = registered + [f"FILLER{i}" for i in range(n_fillers)]
    ids = [f"ENSG{index:011d}" for index in range(1, len(symbols) + 1)]
    rng = np.random.default_rng(seed)
    base = rng.uniform(1.0, 10.0, size=len(symbols))
    senmayo = np.array([symbol in set(SENMAYO_HUMAN) for symbol in symbols])
    soft = hv.parse_geo_samples(soft_text())
    members, member_info = {}, {}
    for sample in soft["samples"]:
        factor = np.ones(len(symbols))
        if sample["condition"] != hv.BASELINE:
            factor[senmayo] = shift
        counts = rng.poisson(np.exp2(base) * factor * rng.uniform(0.9, 1.1))
        lines = [f"{identifier}\t{int(count)}" for identifier, count in zip(ids, counts)]
        lines += ["__no_feature\t100", "__ambiguous\t200", "__too_low_aQual\t0", "__not_aligned\t0",
                  "__alignment_not_unique\t5000"]
        raw = gzip.compress(("\n".join(lines) + "\n").encode(), mtime=0)
        members[sample["count_file"]] = raw
        member_info[sample["count_file"]] = {"sha256": sha(raw), "bytes": len(raw)}
    raw_dir, snapshot, work = Path(directory) / "raw", Path(directory) / "snapshot", Path(directory)
    raw_dir.mkdir()
    snapshot.mkdir()
    tar_path = raw_dir / "GSE160356_RAW.tar"
    with tarfile.open(tar_path, "w") as archive:
        for name, raw in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(raw)
            archive.addfile(info, io.BytesIO(raw))
    hgnc = work / "hgnc_mapping.tsv.gz"
    with gzip.open(hgnc, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["ensembl_gene_id", "symbol"])
        writer.writerows(zip(ids, symbols))
    hv.write_json(snapshot / "samples.json", soft)
    (snapshot / "PLAN.md").write_text("Synthetic frozen plan.\n", encoding="utf-8")
    manifest = {
        "raw_archive": {"url": "https://example.invalid/GSE160356_RAW.tar", "bytes": tar_path.stat().st_size,
                        "sha256": hv.sha256_file(tar_path)},
        "members": member_info, "hgnc": {"sha256": hv.sha256_file(hgnc)},
        "snapshots": {"samples.json": hv.sha256_file(snapshot / "samples.json"),
                      "PLAN.md": hv.sha256_file(snapshot / "PLAN.md")}}
    hv.write_json(snapshot / "sources.json", manifest)
    return {"raw": raw_dir, "snapshot": snapshot, "hgnc": hgnc, "tar": tar_path, "members": members,
            "samples": soft["samples"], "manifest": manifest, "output": work / "artifacts"}


CLEAN = {"revision": "a" * 40, "tracked_changes": False}


class ParseTests(unittest.TestCase):
    def test_a_well_formed_series_is_parsed_with_absent_clones_left_absent(self):
        parsed = hv.parse_geo_samples(soft_text())
        self.assertEqual(len(parsed["samples"]), 9)
        by_condition = {}
        for sample in parsed["samples"]:
            by_condition.setdefault(sample["condition"], []).append(sample["accession"])
            self.assertIsNone(sample["donor_id"])
            self.assertIsNone(sample["clone_id"])
            self.assertIsNone(sample["experiment_block"])
        self.assertEqual({k: len(v) for k, v in by_condition.items()},
                         {"early_passage": 3, "late_passage": 3, "etoposide": 3})
        self.assertEqual(parsed["series"]["title"], "Synthetic fixture")

    def test_unexpected_metadata_is_refused(self):
        good = [sample_block(prefix, i + 1, accession) for prefix, items in GSMS.items()
                for i, accession in enumerate(items)]
        cases = {
            "title": [sample_block("EP_CT", 1, "GSM9000001", title="EP_CT_9")] + good[1:],
            "treatment": [sample_block("EP_CT", 1, "GSM9000001", treatment="Late Passage")] + good[1:],
            "organism": [sample_block("EP_CT", 1, "GSM9000001", organism="Mus musculus")] + good[1:],
            "foreign count file": [sample_block("EP_CT", 1, "GSM9000001",
                                                member="GSM1234567_x.genename.htcounts.txt.gz")] + good[1:],
            "duplicate accession": good + [good[0]],
        }
        for label, blocks in cases.items():
            with self.subTest(label=label), self.assertRaises(ScoreError):
                hv.parse_geo_samples(soft_text(blocks))
        with self.assertRaises(ScoreError):
            hv.parse_geo_samples("^SERIES = GSE160356\n")


class InputTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.fx = build_fixture(self.temp.name, n_fillers=40)

    def test_exactly_the_pinned_members_are_read(self):
        members = hv.extract_members(self.fx["tar"], self.fx["manifest"]["members"])
        self.assertEqual(members, self.fx["members"])

    def test_unexpected_missing_or_altered_members_are_refused(self):
        expected = self.fx["manifest"]["members"]
        name = next(iter(expected))
        with self.assertRaisesRegex(ScoreError, "hash"):
            hv.extract_members(self.fx["tar"], {**expected, name: {"sha256": "0" * 64, "bytes": 1}})
        trimmed = {key: value for key, value in expected.items() if key != name}
        with self.assertRaisesRegex(ScoreError, "unexpected member"):
            hv.extract_members(self.fx["tar"], trimmed)
        with self.assertRaisesRegex(ScoreError, "missing pinned members"):
            hv.extract_members(self.fx["tar"], {**expected, "GSM0000000_extra.genename.htcounts.txt.gz":
                                                {"sha256": "0" * 64, "bytes": 1}})

    def test_a_non_regular_or_oversize_member_is_refused(self):
        expected = self.fx["manifest"]["members"]
        name = next(iter(expected))
        link_tar = Path(self.temp.name) / "link.tar"
        with tarfile.open(link_tar, "w") as archive:
            info = tarfile.TarInfo(name)
            info.type = tarfile.SYMTYPE
            info.linkname = "/etc/passwd"
            archive.addfile(info)
        with self.assertRaisesRegex(ScoreError, "bounded regular file"):
            hv.extract_members(link_tar, {name: expected[name]})
        with patch.object(hv, "MAX_MEMBER_BYTES", 10):
            with self.assertRaisesRegex(ScoreError, "bounded regular file"):
                hv.extract_members(self.fx["tar"], expected)

    def test_assembly_is_deterministic_and_column_aligned(self):
        once = hv.assemble_counts(self.fx["members"], self.fx["samples"])
        twice = hv.assemble_counts(self.fx["members"], self.fx["samples"])
        self.assertEqual(once, twice)
        lines = gzip.decompress(once).decode().splitlines()
        self.assertEqual(lines[0].split("\t"), ["gene", "ep1", "ep2", "ep3", "eto1", "eto2", "eto3", "lp1", "lp2", "lp3"])
        for sample in self.fx["samples"]:
            column = lines[0].split("\t").index(sample["column"])
            original = gzip.decompress(self.fx["members"][sample["count_file"]]).decode().splitlines()
            self.assertEqual([line.split("\t")[column] for line in lines[1:]], [row.split("\t")[1] for row in original])

    def test_assembly_refuses_misaligned_duplicate_and_ragged_tables(self):
        members = dict(self.fx["members"])
        name = self.fx["samples"][0]["count_file"]
        rows = gzip.decompress(members[name]).decode().splitlines()
        cases = {
            "reordered": "\n".join(reversed(rows)) + "\n",
            "duplicate": "\n".join([rows[0], *rows]) + "\n",
            "ragged": "\n".join([rows[0] + "\t1", *rows[1:]]) + "\n",
            "blank identifier": "\n".join(["\t5", *rows[1:]]) + "\n",
        }
        for label, text in cases.items():
            with self.subTest(label=label), self.assertRaises(ScoreError):
                hv.assemble_counts({**members, name: gzip.compress(text.encode(), mtime=0)}, self.fx["samples"])
        with self.assertRaises(ScoreError):
            hv.assemble_counts({**members, name: b"not gzip"}, self.fx["samples"])


class ProtocolTests(unittest.TestCase):
    def test_a_baseline_library_never_helps_choose_its_own_controls(self):
        seen = []

        def fit(rows):
            seen.append(tuple(rows))
            return rows

        def transform(fit_rows):
            return np.array([float(len(fit_rows))] * 9) + np.arange(9) * 0.0

        primary, sensitivity, full = hv.two_protocol_scores([0, 1, 2], fit, transform)
        self.assertEqual(seen, [(0, 1, 2), (1, 2), (0, 2), (0, 1)])
        self.assertEqual(list(sensitivity), [3.0] * 9)
        self.assertEqual(list(primary), [2.0, 2.0, 2.0] + [3.0] * 6)
        self.assertEqual(full, [0, 1, 2])

    def test_primary_and_sensitivity_agree_on_every_held_out_library(self):
        with tempfile.TemporaryDirectory() as directory:
            fx = build_fixture(directory, shift=1.5, n_fillers=600)
            members = fx["members"]
            counts = hv.assemble_counts(members, fx["samples"])
            path = Path(directory) / "counts.tsv.gz"
            path.write_bytes(counts)
            from senescore.public_validation import read_counts
            with gzip.open(fx["hgnc"], "rt") as handle:
                mapping = {row["ensembl_gene_id"]: row["symbol"] for row in csv.DictReader(handle, delimiter="\t")}
            matrix, _raw, names, ordered, _stats = read_counts(path, mapping, fx["samples"])
        baseline = [i for i, s in enumerate(ordered) if s["condition"] == hv.BASELINE]
        key = "senmayo"
        primary, sensitivity, _fit = hv.two_protocol_scores(
            baseline, lambda rows: hv._fit_for(key, matrix, names, ordered, rows, "1" * 64),
            lambda fit: hv._transform(key, matrix, names, fit))
        held_out = [i for i in range(len(ordered)) if i not in baseline]
        np.testing.assert_allclose(primary[held_out], sensitivity[held_out])
        self.assertFalse(np.allclose(primary[baseline], sensitivity[baseline]))

    def test_held_out_expression_cannot_change_fitted_controls(self):
        rng = np.random.default_rng(1)
        names = [f"G{i}" for i in range(400)]
        matrix = rng.normal(5, 2, size=(9, 400))
        signature = names[:30]
        before = fit_module_score(matrix, names, signature, [0, 1, 2], seed=0, n_ctrl=5, n_bins=20,
                                  gene_set_source_status="custom_unverified", citation="synthetic")
        altered = matrix.copy()
        altered[3:] += 50.0
        after = fit_module_score(altered, names, signature, [0, 1, 2], seed=0, n_ctrl=5, n_bins=20,
                                 gene_set_source_status="custom_unverified", citation="synthetic")
        self.assertEqual(before["control_genes"], after["control_genes"])

    def test_empirical_rank_includes_the_target_and_is_two_sided(self):
        self.assertAlmostEqual(hv.empirical_rank(1.0, [0.1, -0.2, 0.3]), 1 / 4)
        self.assertAlmostEqual(hv.empirical_rank(-1.0, [0.1, -1.0, 1.0]), 3 / 4)
        self.assertAlmostEqual(hv.empirical_rank(0.0, [0.1, 0.2]), 1.0)

    def test_the_verdict_rule_is_exactly_the_prespecified_one(self):
        def contrasts(late, etoposide):
            return {"late_passage_vs_early_passage": {"directional_auroc": late},
                    "etoposide_vs_early_passage": {"directional_auroc": etoposide}}

        low = {"late_passage_vs_early_passage": 0.0099, "etoposide_vs_early_passage": 0.05}
        high = {"late_passage_vs_early_passage": 0.0099, "etoposide_vs_early_passage": 0.0501}
        self.assertEqual(hv.verdict(contrasts(1.0, 1.0), low, False), "pass")
        self.assertEqual(hv.verdict(contrasts(1.0, 1.0), high, False),
                         "ranking_only_not_distinguishable_from_random_panels")
        for late, etoposide in ((1.0, 0.9444), (0.5, 1.0), (0.0, 0.0), (0.9999, 1.0)):
            self.assertEqual(hv.verdict(contrasts(late, etoposide), low, False), "fail")
        self.assertEqual(hv.verdict(contrasts(1.0, 1.0), low, True), "not_scored")
        self.assertEqual(hv.verdict(None, low, False), "not_scored")


class EndToEndTests(unittest.TestCase):
    NOW = dt.datetime(2026, 10, 8, 12, 0, tzinfo=dt.timezone.utc)

    @classmethod
    def setUpClass(cls):
        # One full-scale run (all 100 prespecified random panels) is shared by the tests that only
        # read its outputs; the others use a smaller random-panel count so the suite stays quick.
        cls.temp = tempfile.TemporaryDirectory()
        cls.fx = build_fixture(cls.temp.name, shift=3.0)
        with patch.object(hv, "git_state", return_value=CLEAN):
            cls.results = hv.evaluate(cls.fx["raw"], cls.fx["snapshot"], cls.fx["output"], cls.fx["hgnc"],
                                      now=cls.NOW)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def run_small(self, seeds=range(1000, 1030), **fixture):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        fx = build_fixture(temp.name, **fixture)
        with patch.object(hv, "git_state", return_value=CLEAN), patch.object(hv, "RANDOM_SEEDS", seeds):
            results = hv.evaluate(fx["raw"], fx["snapshot"], fx["output"], fx["hgnc"], now=self.NOW)
        return fx, results

    def test_a_strong_synthetic_shift_passes_and_a_null_does_not(self):
        strong = self.results
        self.assertEqual(strong["verdict"]["value"], "pass", strong["verdict"])
        for name in hv.PRIMARY_CONTRASTS:
            self.assertEqual(strong["verdict"]["senmayo_primary_directional_auroc"][name], 1.0)
            self.assertLessEqual(strong["verdict"]["senmayo_empirical_ranks"][name], hv.RANK_THRESHOLD)
        _fx, null = self.run_small(shift=1.0)
        self.assertNotEqual(null["verdict"]["value"], "pass")

    def test_outputs_label_the_split_and_account_for_leakage(self):
        fx, results = self.fx, self.results
        self.assertEqual(results["split"]["baseline_libraries"], sorted(GSMS["EP_CT"]))
        self.assertEqual(sorted(results["split"]["held_out_libraries"]), sorted(GSMS["LP_CT"] + GSMS["ETO"]))
        self.assertEqual([fold["scored"] for fold in results["split"]["leave_one_out_baseline_folds"]],
                         sorted(GSMS["EP_CT"]))
        for fold in results["split"]["leave_one_out_baseline_folds"]:
            self.assertNotIn(fold["scored"], fold["train"])
            self.assertEqual(len(fold["train"]), 2)
        self.assertEqual(results["random_controls"]["n"], 100)
        self.assertEqual(results["code_revision"], "a" * 40)
        self.assertEqual(results["plan_sha256"], fx["manifest"]["snapshots"]["PLAN.md"])
        receipt = json.loads((fx["snapshot"] / "evaluation_receipt.json").read_text())
        leakage = {item["id"]: item["count"] for item in receipt["reported_leakage"]}
        self.assertEqual(leakage, {"held_out_libraries_used_to_fit_controls": 0,
                                   "registered_panel_genes_in_random_panels": 0})
        for name in ("results.json", "random_scores.csv", "evaluation_receipt.json"):
            self.assertTrue((fx["snapshot"] / name).is_file(), name)
        self.assertTrue((fx["output"] / "senmayo_fit.json").is_file())

    def test_every_registered_panel_is_scored_under_both_protocols(self):
        self.assertEqual(set(self.results["panels"]), {*GENE_SETS, "fridman_signed"})
        self.assertEqual(self.results["refusals"], {})
        for panel in self.results["panels"].values():
            for protocol in ("primary", "sensitivity_all_baseline_fit"):
                self.assertEqual(set(panel[protocol]["contrasts"]),
                                 {"late_passage_vs_early_passage", "etoposide_vs_early_passage",
                                  "etoposide_vs_late_passage"})
                self.assertEqual(len(panel[protocol]["scores_by_accession"]), 9)
                for contrast in panel[protocol]["contrasts"].values():
                    self.assertEqual(contrast["permutation_assignments"], 20)

    def test_the_receipt_follows_the_generic_contract(self):
        receipt = json.loads((self.fx["snapshot"] / "evaluation_receipt.json").read_text())
        self.assertEqual(receipt["schema"], "regen-workbench/evaluation-receipt/1")
        self.assertEqual(set(receipt["producer"]), {"tool", "version", "code_revision"})
        self.assertRegex(receipt["producer"]["code_revision"], r"^[0-9a-f]{40}$")
        self.assertTrue(all(len(h) == 64 for h in receipt["input_sha256"]))
        self.assertEqual(len(receipt["input_sha256"]), 10)  # nine count files and the HGNC snapshot
        self.assertEqual(receipt["folds"], [{"train_groups": sorted(GSMS["EP_CT"]),
                                             "test_groups": sorted(GSMS["LP_CT"] + GSMS["ETO"])}])
        self.assertIn(hv.PRIMARY_METRIC, receipt["metric_names"])
        self.assertIn(hv.RANDOM_BASELINE, receipt["model_names"])
        self.assertIsNotNone(dt.datetime.fromisoformat(receipt["created_utc"]).utcoffset())
        self.assertIn("leave-one-out", receipt["baseline_scoring"]["protocol"])

    @unittest.skipUnless((ROOT.parent / "regen-workbench" / "tools" / "frozen_evaluation.py").is_file(),
                         "sibling regen-workbench checkout not present")
    def test_the_receipt_is_accepted_by_the_workbench_adapter_when_it_is_present(self):
        sys.path.insert(0, str(ROOT.parent / "regen-workbench" / "tools"))
        try:
            import frozen_evaluation as fe
        finally:
            sys.path.pop(0)
        normalized = fe.adapt_evaluation_receipt(
            json.loads((self.fx["snapshot"] / "evaluation_receipt.json").read_text()))
        self.assertEqual(normalized["folds"][0]["train_groups"], sorted(GSMS["EP_CT"]))

    def test_it_refuses_a_changed_plan_a_dirty_tree_and_tampered_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            fx = build_fixture(directory, n_fillers=40)
            (fx["snapshot"] / "PLAN.md").write_text("Edited after the freeze.\n", encoding="utf-8")
            with patch.object(hv, "git_state", return_value=CLEAN), self.assertRaisesRegex(ScoreError, "PLAN.md"):
                hv.evaluate(fx["raw"], fx["snapshot"], fx["output"], fx["hgnc"])
        with tempfile.TemporaryDirectory() as directory:
            fx = build_fixture(directory, n_fillers=700)
            dirty = {"revision": "b" * 40, "tracked_changes": True}
            with patch.object(hv, "git_state", return_value=dirty), self.assertRaisesRegex(ScoreError, "differ"):
                hv.evaluate(fx["raw"], fx["snapshot"], fx["output"], fx["hgnc"])
            self.assertFalse((fx["snapshot"] / "results.json").exists())
        with tempfile.TemporaryDirectory() as directory:
            fx = build_fixture(directory, n_fillers=40)
            fx["tar"].write_bytes(fx["tar"].read_bytes() + b"x")
            with self.assertRaisesRegex(ScoreError, "hash mismatch"):
                hv.evaluate(fx["raw"], fx["snapshot"], fx["output"], fx["hgnc"])

    def test_validate_only_reports_structure_and_calculates_no_score(self):
        with tempfile.TemporaryDirectory() as directory:
            fx = build_fixture(directory, n_fillers=60)
            with patch.object(hv, "git_state", side_effect=AssertionError("must not look at git")):
                report = hv.evaluate(fx["raw"], fx["snapshot"], fx["output"], fx["hgnc"], validate_only=True)
            self.assertEqual(set(report), {"accession", "validated", "libraries", "ensembl_rows", "mapped_rows",
                                           "mapped_symbols", "combined_counts_sha256"})
            self.assertEqual(report["libraries"], 9)
            self.assertFalse((fx["snapshot"] / "results.json").exists())
            self.assertFalse(fx["output"].exists())

    def test_a_refused_panel_is_a_recorded_result_not_a_crash(self):
        # Only a few registered genes appear in the table, so every registered panel falls below the
        # 60% coverage requirement and is recorded as a refusal.
        with tempfile.TemporaryDirectory() as directory:
            fx = build_fixture(directory, n_fillers=1500)
            sparse = Path(directory) / "sparse.tsv.gz"
            with gzip.open(fx["hgnc"], "rt") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            keep = {row["symbol"] for row in rows[:5]} | {f"FILLER{i}" for i in range(1500)}
            with gzip.open(sparse, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["ensembl_gene_id", "symbol"])
                writer.writerows((row["ensembl_gene_id"], row["symbol"] if row["symbol"] in keep else "UNRELATED"
                                  ) for row in rows)
            manifest = json.loads((fx["snapshot"] / "sources.json").read_text())
            manifest["hgnc"]["sha256"] = hv.sha256_file(sparse)
            hv.write_json(fx["snapshot"] / "sources.json", manifest)
            with patch.object(hv, "git_state", return_value=CLEAN), patch.object(hv, "RANDOM_SEEDS", range(1000, 1010)):
                results = hv.evaluate(fx["raw"], fx["snapshot"], fx["output"], sparse, now=self.NOW)
        self.assertEqual(results["verdict"]["value"], "not_scored")
        self.assertIn("senmayo", results["refusals"])


SNAPSHOT = ROOT / "validation" / "GSE160356"


@unittest.skipUnless((SNAPSHOT / "results.json").is_file(), "no committed results to check")
class CommittedResultTests(unittest.TestCase):
    """Integrity checks on the committed real-data results, with no download and no raw counts."""

    @classmethod
    def setUpClass(cls):
        cls.sources = json.loads((SNAPSHOT / "sources.json").read_text(encoding="utf-8"))
        cls.results = json.loads((SNAPSHOT / "results.json").read_text(encoding="utf-8"))
        cls.receipt = json.loads((SNAPSHOT / "evaluation_receipt.json").read_text(encoding="utf-8"))

    def test_the_plan_and_sample_snapshots_are_the_frozen_ones(self):
        for name in ("PLAN.md", "samples.json"):
            self.assertEqual(hv.sha256_file(SNAPSHOT / name), self.sources["snapshots"][name], name)
        self.assertEqual(self.results["plan_sha256"], self.sources["snapshots"]["PLAN.md"])

    def test_the_run_came_from_one_clean_committed_revision(self):
        self.assertRegex(self.results["code_revision"], r"^[0-9a-f]{40}$")
        self.assertIs(self.results["tracked_tree_changes"], False)
        self.assertEqual(self.receipt["producer"]["code_revision"], self.results["code_revision"])
        self.assertEqual(self.receipt["created_utc"], self.results["created_utc"])

    def test_the_stored_verdict_follows_from_the_stored_numbers(self):
        senmayo = self.results["panels"]["senmayo"]["primary"]["contrasts"]
        ranks = self.results["verdict"]["senmayo_empirical_ranks"]
        self.assertEqual(self.results["verdict"]["value"], hv.verdict(senmayo, ranks, "senmayo" in self.results["refusals"]))

    def test_the_stored_ranks_recompute_from_the_committed_random_scores(self):
        with (SNAPSHOT / "random_scores.csv").open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual([int(row["seed"]) for row in rows], list(range(1000, 1100)))
        for name, _positive, _negative in hv.CONTRASTS:
            deltas = [float(row[f"{name}__mean_difference"]) for row in rows]
            target = self.results["panels"]["senmayo"]["primary"]["contrasts"][name]["mean_difference"]
            self.assertAlmostEqual(hv.empirical_rank(target, deltas),
                                   self.results["verdict"]["senmayo_empirical_ranks"][name])

    def test_every_pinned_input_hash_is_in_the_receipt_and_nothing_else_is(self):
        pinned = {info["sha256"] for info in self.sources["members"].values()} | {self.sources["hgnc"]["sha256"]}
        self.assertEqual(set(self.receipt["input_sha256"]), pinned)

    def test_no_held_out_library_was_allowed_to_influence_a_fit(self):
        leakage = {item["id"]: item["count"] for item in self.receipt["reported_leakage"]}
        self.assertEqual(leakage["held_out_libraries_used_to_fit_controls"], 0)
        baseline = set(self.results["split"]["baseline_libraries"])
        for fold in self.results["split"]["leave_one_out_baseline_folds"]:
            self.assertTrue(set(fold["train"]) < baseline)
            self.assertNotIn(fold["scored"], fold["train"])


if __name__ == "__main__":
    unittest.main()

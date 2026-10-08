"""Checks for the SenMayo gene-removal robustness run, on a SYNTHETIC table.

No real GSE160356 score is calculated here. The fixture has the real layout (nine libraries,
three per condition, all registered SenMayo symbols plus filler genes) with a made-up signal.
"""

import csv
import datetime as dt
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from senescore import hrmec_robustness as hr
from senescore.genes import GENE_SETS, ORTHOGONAL

SENMAYO = list(GENE_SETS["senmayo"]["symbols"])
CONDITIONS = ["early_passage"] * 3 + ["late_passage"] * 3 + ["etoposide"] * 3
ACCESSIONS = [f"GSM90000{i:02d}" for i in range(9)]


def write_fixture(root: Path):
    rng = np.random.default_rng(5)
    names = SENMAYO + list(ORTHOGONAL) + [f"FILL{i}" for i in range(400)]
    base = rng.uniform(1, 9, size=len(names))
    matrix = base + rng.normal(0, 0.15, size=(9, len(names)))
    induced = [i for i, c in enumerate(CONDITIONS) if c != "early_passage"]
    matrix[np.ix_(induced, range(len(SENMAYO)))] += 1.0  # planted shift on every panel gene
    table = root / "table.csv"
    with table.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample_id", *names])
        for accession, row in zip(ACCESSIONS, matrix):
            writer.writerow([accession, *row])
    samples = root / "samples.json"
    samples.write_text(json.dumps({"samples": [{"accession": a, "condition": c}
                                               for a, c in zip(ACCESSIONS, CONDITIONS)]}))
    groups = {}
    for i, c in enumerate(CONDITIONS):
        groups.setdefault(c, []).append(i)
    accessions, loaded_names, loaded = hr.load_table(table)
    primary, _ = hr.score_panel(loaded, loaded_names, accessions, groups, SENMAYO, registered=True)
    contrasts = hr.contrasts_for(primary, groups, permutation=False)
    results = root / "results.json"
    results.write_text(json.dumps({"panels": {"senmayo": {"primary": {
        "scores_by_accession": {a: float(v) for a, v in zip(ACCESSIONS, primary)},
        "contrasts": {n: {"mean_difference": contrasts[n]["mean_difference"]} for n in hr.PRIMARY_CONTRASTS}}}}}))
    random_scores = root / "random_scores.csv"
    with random_scores.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([f"{n}__mean_difference" for n, *_ in hr.CONTRASTS])
        for value in np.linspace(-0.2, 0.2, 100):
            writer.writerow([value] * len(hr.CONTRASTS))
    return table, samples, results, random_scores, groups, (accessions, loaded_names, loaded)


class RemovalTests(unittest.TestCase):
    def test_a_draw_is_deterministic_sorted_distinct_and_from_the_present_genes(self):
        present = SENMAYO[:124]
        first = hr.removal_draw(present, 2000)
        self.assertEqual(first, hr.removal_draw(list(reversed(present)), 2000))
        self.assertEqual(len(first), hr.N_REMOVED)
        self.assertEqual(first, sorted(first))
        self.assertEqual(len(set(first)), hr.N_REMOVED)
        self.assertTrue(set(first) <= set(present))
        self.assertNotEqual(first, hr.removal_draw(present, 2001))

    def test_a_draw_cannot_empty_the_panel(self):
        with self.assertRaises(ValueError):
            hr.removal_draw(["A", "B", "C"], 1, n_removed=3)
        with self.assertRaises(ValueError):
            hr.removal_draw(["A", "B", "C"], 1, n_removed=0)

    def test_the_frozen_seed_range_and_removal_size(self):
        self.assertEqual((DRAWS := hr.DRAW_SEEDS).start, 2000)
        self.assertEqual(len(DRAWS), 1000)
        self.assertEqual(hr.N_REMOVED, 10)


class RunTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        (self.root / "out").mkdir()
        (self.root / "out" / "PLAN.md").write_text("plan\n")
        (self.table, self.samples, self.results, self.random,
         self.groups, self.loaded) = write_fixture(self.root)

    def test_removed_genes_stay_out_of_the_control_pool_and_the_panel(self):
        accessions, names, matrix = self.loaded
        removed = hr.removal_draw(SENMAYO, 2000)
        kept = [g for g in SENMAYO if g not in removed]
        baseline = self.groups["early_passage"]
        _, coverage = hr.score_panel(matrix, names, accessions, self.groups, kept, removed)
        self.assertEqual(coverage, 1.0)
        fit = hr.fit_module_score(matrix, names, kept, baseline, seed=0, n_ctrl=5, n_bins=20,
                                  block_from_controls=set(ORTHOGONAL) | set(removed),
                                  orthogonal=ORTHOGONAL, gene_set_source_status="custom_unverified")
        controls = {c for chosen in fit["control_genes"].values() for c in chosen}
        self.assertFalse(controls & set(removed))
        self.assertFalse(controls & set(ORTHOGONAL))
        self.assertEqual(set(fit["control_genes"]), set(kept))

    def test_a_small_run_writes_results_and_draws_and_a_strong_planted_shift_keeps_pass(self):
        results = hr.run(self.table, self.samples, self.results, self.random, self.root / "out",
                         draw_seeds=range(2000, 2004), allow_dirty=True,
                         now=dt.datetime(2026, 10, 8, tzinfo=dt.timezone.utc))
        summary = results["summary"]
        self.assertEqual(summary["n_draws"], 4)
        self.assertEqual(summary["draws_keeping_pass"] + summary["draws_losing_pass"], 4)
        self.assertEqual(summary["draws_keeping_pass"], 4)
        self.assertEqual(results["zero_removal_control_max_abs_gap"], 0.0)
        self.assertEqual(results["present_senmayo_genes"], len(SENMAYO))
        with (self.root / "out" / "draws.csv").open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual([r["seed"] for r in rows], ["2000", "2001", "2002", "2003"])
        self.assertEqual(len(rows[0]["removed_genes"].split()), hr.N_REMOVED)

    def test_the_control_must_reproduce_the_committed_scores(self):
        data = json.loads(self.results.read_text())
        first = ACCESSIONS[0]
        data["panels"]["senmayo"]["primary"]["scores_by_accession"][first] += 0.01
        self.results.write_text(json.dumps(data))
        with self.assertRaisesRegex(hr.ScoreError, "zero-removal control"):
            hr.run(self.table, self.samples, self.results, self.random, self.root / "out",
                   draw_seeds=range(2000, 2001), allow_dirty=True)
        self.assertFalse((self.root / "out" / "results.json").exists())

    def test_a_dirty_tree_is_refused_for_the_frozen_run(self):
        from unittest.mock import patch
        with patch.object(hr, "git_state", return_value={"revision": "0" * 40, "tracked_changes": ["x"]}):
            with self.assertRaisesRegex(hr.ScoreError, "commit before the frozen run"):
                hr.run(self.table, self.samples, self.results, self.random, self.root / "out",
                       draw_seeds=range(2000, 2001))

    def test_a_losing_verdict_is_counted_not_hidden(self):
        # Random-panel deltas wider than any real shift put every rank above the threshold.
        with self.random.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow([f"{n}__mean_difference" for n, *_ in hr.CONTRASTS])
            for _ in range(100):
                writer.writerow([50.0] * len(hr.CONTRASTS))
        results = hr.run(self.table, self.samples, self.results, self.random, self.root / "out",
                         draw_seeds=range(2000, 2002), allow_dirty=True)
        self.assertEqual(results["summary"]["draws_losing_pass"], 2)
        self.assertIn("ranking_only_not_distinguishable_from_random_panels", results["summary"]["verdict_counts"])


if __name__ == "__main__":
    unittest.main()

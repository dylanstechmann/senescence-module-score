"""Random removal of SenMayo genes from the GSE160356 check (frozen plan: validation/GSE160356/robustness/PLAN.md).

Reads the transformed log2(CPM+1) table the earlier run wrote, removes 10 randomly chosen
SenMayo genes per draw (seeds 2000..2999), scores the nine libraries under the earlier primary
protocol and applies the earlier verdict rule unchanged. A zero-removal control must reproduce the
committed SenMayo scores before any draw is trusted.

This is a robustness check on the same nine libraries, not a new test on independent data.

    python -m senescore.hrmec_robustness --out validation/GSE160356/robustness
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path

import numpy as np

from senescore.genes import GENE_SETS, ORTHOGONAL
from senescore.hrmec_validation import (
    BASELINE, CONTRASTS, PRIMARY_CONTRASTS, contrasts_for, empirical_rank, git_state,
    two_protocol_scores, verdict,
)
from senescore.score import ScoreError, fit_module_score, transform_module_score

DRAW_SEEDS = range(2000, 3000)
N_REMOVED = 10
CONTROL_TOLERANCE = 1e-9
HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]
VALIDATION = ROOT / "validation" / "GSE160356"


def load_table(path: Path):
    """Read ``sample_id, gene...`` rows. Returns accessions, gene names and the matrix."""
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        rows = list(reader)
    accessions = [row[0] for row in rows]
    matrix = np.array([[float(value) for value in row[1:]] for row in rows], dtype=float)
    if not np.all(np.isfinite(matrix)):
        raise ScoreError("expression table has non-finite values")
    return accessions, header[1:], matrix


def removal_draw(present: list[str], seed: int, n_removed: int = N_REMOVED) -> list[str]:
    """The genes removed by one draw: uniform, without replacement, from the sorted present genes."""
    if not 0 < n_removed < len(present):
        raise ValueError("n_removed must leave a non-empty panel")
    rng = np.random.default_rng(seed)
    return sorted(rng.choice(sorted(present), size=n_removed, replace=False).tolist())


def score_panel(matrix, names, accessions, groups, symbols, removed=(), *, registered=False):
    """Primary-protocol scores for one panel; removed genes also stay out of the control pool."""
    baseline = groups[BASELINE]
    entry = GENE_SETS["senmayo"]
    blocked = set(ORTHOGONAL) | set(removed)

    def fit(rows):
        if registered:
            return fit_module_score(matrix, names, entry["symbols"], rows, seed=0, n_ctrl=5, n_bins=20,
                                    block_from_controls=blocked, orthogonal=entry["orthogonal"],
                                    citation=entry["citation"], gene_set_key="senmayo",
                                    gene_set_source_status=entry["source_status"])
        return fit_module_score(matrix, names, list(symbols), rows, seed=0, n_ctrl=5, n_bins=20,
                                block_from_controls=blocked, orthogonal=entry["orthogonal"],
                                citation=entry["citation"], gene_set_source_status="custom_unverified")

    primary, _sensitivity, full_fit = two_protocol_scores(
        baseline, fit, lambda f: np.asarray(transform_module_score(matrix, names, f)["scores"], dtype=float))
    return primary, full_fit["coverage"]


def sha256_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(table_path: Path, samples_path: Path, results_path: Path, random_scores_path: Path,
        out_dir: Path, *, draw_seeds=DRAW_SEEDS, allow_dirty=False, now=None) -> dict:
    state = git_state(ROOT)
    if state["tracked_changes"] and not allow_dirty:
        raise ScoreError("tracked files differ from the committed revision; commit before the frozen run")
    accessions, names, matrix = load_table(table_path)
    condition = {s["accession"]: s["condition"] for s in json.loads(Path(samples_path).read_text("utf-8"))["samples"]}
    groups: dict[str, list[int]] = {}
    for index, accession in enumerate(accessions):
        groups.setdefault(condition[accession], []).append(index)
    committed = json.loads(Path(results_path).read_text("utf-8"))
    reference = committed["panels"]["senmayo"]["primary"]["scores_by_accession"]

    full, _coverage = score_panel(matrix, names, accessions, groups, GENE_SETS["senmayo"]["symbols"],
                                  registered=True)
    control_gap = max(abs(float(full[i]) - reference[a]) for i, a in enumerate(accessions))
    if control_gap > CONTROL_TOLERANCE:
        raise ScoreError(f"zero-removal control differs from the committed scores by {control_gap:.3g}")

    with Path(random_scores_path).open(encoding="utf-8", newline="") as handle:
        random_rows = list(csv.DictReader(handle))
    random_deltas = {name: [float(r[f"{name}__mean_difference"]) for r in random_rows] for name, *_ in CONTRASTS}
    present = [g for g in GENE_SETS["senmayo"]["symbols"] if g in set(names)]

    draws = []
    for seed in draw_seeds:
        removed = removal_draw(present, seed)
        kept = [g for g in present if g not in set(removed)]
        scores, coverage = score_panel(matrix, names, accessions, groups, kept, removed)
        contrasts = contrasts_for(scores, groups, permutation=False)
        ranks = {name: empirical_rank(contrasts[name]["mean_difference"], random_deltas[name])
                 for name in PRIMARY_CONTRASTS}
        draws.append({"seed": seed, "removed": removed, "coverage": coverage,
                      "contrasts": {name: {"mean_difference": contrasts[name]["mean_difference"],
                                           "directional_auroc": contrasts[name]["directional_auroc"]}
                                    for name in PRIMARY_CONTRASTS},
                      "ranks": ranks, "verdict": verdict(contrasts, ranks, False)})

    def quantiles(values):
        return np.quantile(values, [0, .025, .5, .975, 1]).tolist()

    kept_pass = sum(d["verdict"] == "pass" for d in draws)
    summary = {
        "n_draws": len(draws), "n_removed": N_REMOVED, "panel_size_after_removal": len(present) - N_REMOVED,
        "draws_keeping_pass": kept_pass, "draws_losing_pass": len(draws) - kept_pass,
        "draws_with_both_primary_auroc_one": sum(
            all(d["contrasts"][n]["directional_auroc"] == 1.0 for n in PRIMARY_CONTRASTS) for d in draws),
        "verdict_counts": {v: sum(d["verdict"] == v for d in draws) for v in sorted({d["verdict"] for d in draws})},
        "mean_difference_quantiles_min_2.5_median_97.5_max": {
            n: quantiles([d["contrasts"][n]["mean_difference"] for d in draws]) for n in PRIMARY_CONTRASTS},
        "committed_full_panel_mean_difference": {
            n: committed["panels"]["senmayo"]["primary"]["contrasts"][n]["mean_difference"] for n in PRIMARY_CONTRASTS},
        "max_rank_observed": {n: max(d["ranks"][n] for d in draws) for n in PRIMARY_CONTRASTS},
        "smallest_late_passage_shift_draws": [
            {"seed": d["seed"], "removed": d["removed"],
             "mean_difference": d["contrasts"]["late_passage_vs_early_passage"]["mean_difference"],
             "verdict": d["verdict"]}
            for d in sorted(draws, key=lambda d: d["contrasts"]["late_passage_vs_early_passage"]["mean_difference"])[:5]],
    }
    created = (now or dt.datetime.now(dt.timezone.utc)).isoformat(timespec="seconds")
    results = {
        "plan": "validation/GSE160356/robustness/PLAN.md", "plan_sha256": sha256_file(out_dir / "PLAN.md"),
        "created_utc": created, "code_revision": state["revision"], "tracked_tree_changes": state["tracked_changes"],
        "input_table_sha256": sha256_file(table_path), "committed_results_sha256": sha256_file(results_path),
        "random_scores_sha256": sha256_file(random_scores_path),
        "zero_removal_control_max_abs_gap": control_gap, "zero_removal_control_tolerance": CONTROL_TOLERANCE,
        "present_senmayo_genes": len(present), "summary": summary,
        "limits": ["same nine libraries as the earlier check; not independent data",
                   "the 100 random panels have 125 genes and the reduced panels 114, so ranks are approximate",
                   "random removal does not identify which genes drive the shift",
                   "no tissue generalization, rejuvenation or diagnostic claim"],
    }
    (out_dir / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True, allow_nan=False) + "\n",
                                          encoding="utf-8")
    with (out_dir / "draws.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["seed", "removed_genes", "late_minus_early", "late_auroc", "late_rank",
                         "eto_minus_early", "eto_auroc", "eto_rank", "verdict"])
        for d in draws:
            late, eto = (d["contrasts"][n] for n in PRIMARY_CONTRASTS)
            writer.writerow([d["seed"], " ".join(d["removed"]), late["mean_difference"], late["directional_auroc"],
                             d["ranks"][PRIMARY_CONTRASTS[0]], eto["mean_difference"], eto["directional_auroc"],
                             d["ranks"][PRIMARY_CONTRASTS[1]], d["verdict"]])
    return results


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default=str(VALIDATION / "robustness"))
    parser.add_argument("--table", default=str(ROOT / "artifacts" / "GSE160356" / "expression_logCPM.csv"))
    parser.add_argument("--allow-dirty", action="store_true", help="for tests only; the frozen run must be clean")
    args = parser.parse_args(argv)
    out = Path(args.out)
    if (out / "results.json").exists():
        raise SystemExit("results.json exists; the frozen run happens once")
    results = run(Path(args.table), VALIDATION / "samples.json", VALIDATION / "results.json",
                  VALIDATION / "random_scores.csv", out, allow_dirty=args.allow_dirty)
    print(json.dumps(results["summary"], indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

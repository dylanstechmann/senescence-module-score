"""Four fixed gene removals on GSE160356 (frozen plan: validation/GSE160356/targeted/PLAN.md).

    python -m senescore.hrmec_targeted
"""
from __future__ import annotations

import csv
import datetime as dt
import json
from pathlib import Path

from senescore.genes import GENE_SETS
from senescore.hrmec_robustness import ROOT, VALIDATION, load_table, score_panel, sha256_file
from senescore.hrmec_validation import PRIMARY_CONTRASTS, CONTRASTS, contrasts_for, empirical_rank, git_state, verdict
from senescore.score import ScoreError

RECURRING3 = ("PLAT", "CSF2RB", "IL6")
LATE_TOP10 = ("PLAT", "IGFBP5", "SERPINE2", "MMP1", "IGFBP6", "IGFBP3", "IL18", "FAS", "ANGPTL4", "CSF2RB")
ETOPOSIDE_TOP10 = ("MMP1", "SERPINE2", "IGFBP3", "MMP10", "FAS", "ANGPT1", "PLAT", "ANGPTL4", "CXCL8", "CSF2RB")
REMOVALS = {"recurring3": RECURRING3, "late_top10": LATE_TOP10, "etoposide_top10": ETOPOSIDE_TOP10,
            "union_top": tuple(sorted(set(RECURRING3) | set(LATE_TOP10) | set(ETOPOSIDE_TOP10)))}


def run(out_dir: Path, table_path: Path, *, allow_dirty=False, now=None) -> dict:
    state = git_state(ROOT)
    if state["tracked_changes"] and not allow_dirty:
        raise ScoreError("commit before the frozen run")
    accessions, names, matrix = load_table(table_path)
    condition = {s["accession"]: s["condition"]
                 for s in json.loads((VALIDATION / "samples.json").read_text("utf-8"))["samples"]}
    groups: dict[str, list[int]] = {}
    for i, a in enumerate(accessions):
        groups.setdefault(condition[a], []).append(i)
    with (VALIDATION / "random_scores.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    deltas = {n: [float(r[f"{n}__mean_difference"]) for r in rows] for n, *_ in CONTRASTS}
    present = [g for g in GENE_SETS["senmayo"]["symbols"] if g in set(names)]
    results = {}
    for label, removed in REMOVALS.items():
        missing = [g for g in removed if g not in present]
        if missing:
            raise ScoreError(f"{label}: genes not in the present panel: {missing}")
        kept = [g for g in present if g not in removed]
        scores, _ = score_panel(matrix, names, accessions, groups, kept, removed)
        contrasts = contrasts_for(scores, groups, permutation=False)
        ranks = {n: empirical_rank(contrasts[n]["mean_difference"], deltas[n]) for n in PRIMARY_CONTRASTS}
        results[label] = {"removed": list(removed), "panel_size": len(kept),
                          "contrasts": {n: {"mean_difference": contrasts[n]["mean_difference"],
                                            "directional_auroc": contrasts[n]["directional_auroc"]}
                                        for n in PRIMARY_CONTRASTS},
                          "ranks": ranks, "verdict": verdict(contrasts, ranks, False)}
    created = (now or dt.datetime.now(dt.timezone.utc)).isoformat(timespec="seconds")
    payload = {"plan": "validation/GSE160356/targeted/PLAN.md", "plan_sha256": sha256_file(out_dir / "PLAN.md"),
               "created_utc": created, "code_revision": state["revision"],
               "tracked_tree_changes": state["tracked_changes"], "input_table_sha256": sha256_file(table_path),
               "removals_losing_pass": sum(r["verdict"] != "pass" for r in results.values()),
               "removals": results,
               "limits": ["same nine libraries; genes were chosen from these data", "ranks are approximate (smaller panels)"]}
    (out_dir / "results.json").write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
                                          encoding="utf-8")
    return payload


def main():
    out = VALIDATION / "targeted"
    if (out / "results.json").exists():
        raise SystemExit("results.json exists; the frozen run happens once")
    print(json.dumps(run(out, ROOT / "artifacts" / "GSE160356" / "expression_logCPM.csv"), indent=2))


if __name__ == "__main__":
    main()

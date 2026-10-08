"""Post-hoc and descriptive, NOT part of the frozen plan: which removed genes recur in draws that lose `pass`?

Reads draws.csv (written by the frozen run) and prints, for each gene removed at least 40 times,
the share of its draws that lost `pass`. With 1,000 draws and 10 of 124 genes removed each time, a
gene is removed in about 8% of draws; the overall loss rate is the baseline to compare against.
Nothing here changes a verdict, and a gene that recurs is a lead for a separately frozen targeted
removal, not a finding.

    python validation/GSE160356/robustness/posthoc_removed_genes.py
"""
import collections
import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main():
    rows = list(csv.DictReader((HERE / "draws.csv").open(encoding="utf-8", newline="")))
    drawn, lost = collections.Counter(), collections.Counter()
    for row in rows:
        for gene in row["removed_genes"].split():
            drawn[gene] += 1
            lost[gene] += row["verdict"] != "pass"
    overall = sum(row["verdict"] != "pass" for row in rows) / len(rows)
    table = sorted(({"gene": g, "draws_removed": drawn[g], "draws_losing_pass": lost[g],
                     "loss_rate": round(lost[g] / drawn[g], 3)} for g in drawn if drawn[g] >= 40),
                   key=lambda item: (-item["loss_rate"], item["gene"]))
    result = {"note": "post-hoc and descriptive; not part of the frozen plan", "overall_loss_rate": round(overall, 3),
              "top_by_loss_rate": table[:10]}
    (HERE / "posthoc_removed_genes.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

"""Post-hoc, descriptive, NOT part of the prespecified verdict: which SenMayo genes move?

Reads the transformed log2(CPM+1) table written by ``python -m senescore.hrmec_validation``
(``artifacts/GSE160356/expression_logCPM.csv``, ignored by git) and, for each SenMayo gene present,
reports mean(induced) - mean(early passage) in log2(CPM+1). No fit, control or random panel is
involved, nothing here changes a panel or the verdict, and it was written after the verdict existed.

    python validation/GSE160356/posthoc_gene_contributions.py
"""
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from senescore.genes import SENMAYO_HUMAN  # noqa: E402

TABLE = ROOT / "artifacts" / "GSE160356" / "expression_logCPM.csv"
SAMPLES = json.loads((Path(__file__).parent / "samples.json").read_text(encoding="utf-8"))["samples"]
CONDITION = {s["accession"]: s["condition"] for s in SAMPLES}
TITLE = {s["accession"]: s["title"] for s in SAMPLES}


def main():
    with TABLE.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        rows = {row[0]: [float(value) for value in row[1:]] for row in reader}
    genes = header[1:]
    index = {gene: i for i, gene in enumerate(genes)}
    present = [gene for gene in SENMAYO_HUMAN if gene in index]
    baseline = [a for a, c in CONDITION.items() if c == "early_passage"]

    def mean(accessions, gene):
        return sum(rows[a][index[gene]] for a in accessions) / len(accessions)

    result = {"note": "post-hoc and descriptive; log2(CPM+1) difference of group means, no controls", "contrasts": {}}
    for condition in ("late_passage", "etoposide"):
        induced = [a for a, c in CONDITION.items() if c == condition]
        deltas = sorted(((mean(induced, g) - mean(baseline, g), g) for g in present), reverse=True)
        values = [d for d, _ in deltas]
        top = [g for _, g in deltas[:10]]
        remaining = [d for d, g in deltas if g not in top]
        result["contrasts"][f"{condition}_vs_early_passage"] = {
            "n_senmayo_genes_present": len(present),
            "n_genes_higher_in_induced": sum(v > 0 for v in values),
            "n_genes_lower_in_induced": sum(v < 0 for v in values),
            "mean_gene_difference_all": sum(values) / len(values),
            "mean_gene_difference_without_top10": sum(remaining) / len(remaining),
            "top10_higher": [{"gene": g, "difference": round(d, 3)} for d, g in deltas[:10]],
            "top5_lower": [{"gene": g, "difference": round(d, 3)} for d, g in deltas[-5:][::-1]],
        }
    out = Path(__file__).parent / "posthoc_gene_contributions.json"
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

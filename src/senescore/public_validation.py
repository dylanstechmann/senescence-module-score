"""Reproduce the prespecified external LF1 expression check, GEO GSE268487.

This intentionally narrow importer never infers donors, joins separate assays,
changes signature members or fits controls using held-out conditions.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import re
import urllib.request

import numpy as np

from senescore.genes import GENE_SETS, ORTHOGONAL, SENMAYO_HUMAN
from senescore.score import (
    ScoreError, cohen_d, fit_module_score, fit_signed_fridman_score,
    pearson_correlation, random_signature, transform_module_score,
    transform_signed_fridman_score,
)

ACCESSION = "GSE268487"
IMPORTER_VERSION = 3
ENSEMBL = re.compile(r"^ENSG\d{11}(?:\.\d+)?(?:_PAR_Y)?$")
TITLE = re.compile(r"^LF1 cells, (proliferating|quiescent|senescent), replicate ([1-3])$")
CONDITIONS = ("proliferating", "quiescent", "senescent")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")


def parse_geo_samples(text: str) -> list[dict]:
    """Extract exact accession/title/condition facts; absent donors stay absent."""
    records, current = [], None
    for line in text.splitlines():
        if line.startswith("^"):
            if current is not None:
                records.append(current)
            current = {"accession": line.split(" = ", 1)[1], "fields": {}} if line.startswith("^SAMPLE = ") else None
        elif current is not None and line.startswith("!Sample_") and " = " in line:
            key, value = line.split(" = ", 1)
            current["fields"].setdefault(key, []).append(value)
    if current is not None:
        records.append(current)
    samples = []
    for record in records:
        fields = record["fields"]
        titles = fields.get("!Sample_title", [])
        match = TITLE.fullmatch(titles[0]) if len(titles) == 1 else None
        if not match or fields.get("!Sample_organism_ch1") != ["Homo sapiens"]:
            raise ScoreError("unexpected or incomplete GSE268487 sample metadata")
        condition, replicate = match.groups()
        characteristics = fields.get("!Sample_characteristics_ch1", [])
        expected_treatment = {"proliferating": "None", "quiescent": "Contact Inhibition + Serum Starvation",
                              "senescent": "Replicative Exhaustion"}[condition]
        if ("cell line: LF1" not in characteristics or "cell type: lung fibroblasts" not in characteristics
                or f"treatment: {expected_treatment}" not in characteristics):
            raise ScoreError("sample metadata has unexpected LF1 cell type or treatment")
        prefix = {"proliferating": "pro", "quiescent": "qui", "senescent": "sen"}[condition]
        samples.append({
            "accession": record["accession"], "title": titles[0],
            "column": prefix + replicate, "condition": condition,
            "reported_replicate": int(replicate), "cell_line": "LF1",
            "donor_id": None, "experiment_block": None,
            "source_fields": {key: values for key, values in fields.items() if key in (
                "!Sample_organism_ch1", "!Sample_characteristics_ch1",
                "!Sample_treatment_protocol_ch1", "!Sample_extract_protocol_ch1",
                "!Sample_data_processing", "!Sample_relation", "!Sample_status",
            )},
        })
    if not samples:
        raise ScoreError("GEO metadata contains no samples")
    if len({s["column"] for s in samples}) != len(samples) or len({s["accession"] for s in samples}) != len(samples):
        raise ScoreError("duplicate sample accession or condition/replicate")
    return samples


def reduce_hgnc(text: str) -> tuple[dict[str, str], list[str]]:
    """Only approved Ensembl-to-current-symbol mappings, rejecting ambiguity."""
    reader = csv.DictReader(text.splitlines(), delimiter="\t")
    if not {"ensembl_gene_id", "symbol", "status"}.issubset(reader.fieldnames or []):
        raise ScoreError("HGNC snapshot lacks required columns")
    candidates: dict[str, set[str]] = {}
    for row in reader:
        if row["status"] != "Approved" or not row["ensembl_gene_id"]:
            continue
        identifier, symbol = row["ensembl_gene_id"].strip(), row["symbol"].strip().upper()
        if not ENSEMBL.fullmatch(identifier) or not symbol:
            raise ScoreError("invalid HGNC Ensembl/symbol field")
        base = identifier.split(".")[0]
        candidates.setdefault(base, set()).add(symbol)
    ambiguous = sorted(key for key, symbols in candidates.items() if len(symbols) != 1)
    mapping = {key: next(iter(symbols)) for key, symbols in candidates.items() if len(symbols) == 1}
    return mapping, ambiguous


def read_counts(path: Path, mapping: dict[str, str], samples: list[dict]):
    """Return sample-by-symbol logCPM; unmapped ENSG counts remain in denominator."""
    columns_to_samples = {sample["column"]: sample for sample in samples}
    if len(columns_to_samples) != len(samples):
        raise ScoreError("duplicate metadata columns")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            raise ScoreError("count table is empty") from None
        columns = [value.strip() for value in header[1:]]
        if not columns or len(columns) != len(set(columns)) or set(columns) != set(columns_to_samples):
            raise ScoreError("count columns do not match unique GEO sample metadata")
        ordered_samples = [columns_to_samples[column] for column in columns]
        totals = np.zeros(len(columns), dtype=float)
        excluded = np.zeros(len(columns), dtype=float)
        symbol_counts: dict[str, np.ndarray] = {}
        identifiers, ensembl_bases = set(), set()
        statistics = {"ensembl_rows": 0, "mapped_rows": 0, "unmapped_ensembl_rows": 0, "excluded_non_ensembl_rows": 0}
        for row_number, row in enumerate(reader, start=2):
            if len(row) != len(columns) + 1 or not row[0].strip():
                raise ScoreError(f"ragged or blank count row {row_number}")
            identifier = row[0].strip()
            if identifier in identifiers:
                raise ScoreError(f"duplicate count identifier {identifier}")
            identifiers.add(identifier)
            try:
                values = np.array(row[1:], dtype=float)
            except ValueError:
                raise ScoreError(f"nonnumeric count row {row_number}") from None
            if not np.isfinite(values).all() or (values < 0).any() or (values != np.floor(values)).any():
                raise ScoreError(f"counts must be finite nonnegative integers (row {row_number})")
            if not ENSEMBL.fullmatch(identifier):
                if identifier.startswith("ENSG"):
                    raise ScoreError(f"malformed Ensembl identifier {identifier}")
                excluded += values
                statistics["excluded_non_ensembl_rows"] += 1
                continue
            base = identifier.split(".")[0].removesuffix("_PAR_Y")
            locus = base + ("_PAR_Y" if identifier.endswith("_PAR_Y") else "")
            if locus in ensembl_bases:
                raise ScoreError(f"duplicate version-stripped Ensembl identifier {locus}")
            ensembl_bases.add(locus)
            totals += values
            statistics["ensembl_rows"] += 1
            if base not in mapping:
                statistics["unmapped_ensembl_rows"] += 1
                continue
            symbol = mapping[base]
            if not symbol or symbol != symbol.strip().upper():
                raise ScoreError("mapping symbols must be normalized nonblank symbols")
            symbol_counts.setdefault(symbol, np.zeros(len(columns), dtype=float))[:] += values
            statistics["mapped_rows"] += 1
    if not symbol_counts or (totals <= 0).any() or not np.isfinite(totals).all():
        raise ScoreError("count table has no mapped genes or a zero/invalid gene library size")
    names = sorted(symbol_counts)
    raw = np.column_stack([symbol_counts[symbol] for symbol in names])
    matrix = np.log2(raw / totals[:, None] * 1_000_000 + 1)
    statistics.update({"mapped_symbols": len(names), "ensembl_library_counts": totals.tolist(),
                       "excluded_non_ensembl_counts": excluded.tolist()})
    return matrix, raw, names, ordered_samples, statistics


def descriptive_contrast(scores, positive_rows, negative_rows, *, permutation=False):
    positive, negative = np.asarray(scores)[positive_rows], np.asarray(scores)[negative_rows]
    difference = float(positive.mean() - negative.mean())
    comparisons = positive[:, None] - negative[None, :]
    auc = float(((comparisons > 0).sum() + 0.5 * (comparisons == 0).sum()) / comparisons.size)
    result = {"n_positive": len(positive), "n_negative": len(negative),
              "mean_difference": difference, "directional_auroc": auc,
              "cohen_d": None, "cohen_d_note": None}
    try:
        result["cohen_d"] = cohen_d(positive, negative)
    except ScoreError as exc:
        result["cohen_d_note"] = str(exc)
    if permutation:
        pooled = np.concatenate((positive, negative))
        assignments = list(itertools.combinations(range(len(pooled)), len(positive)))
        extreme = 0
        for selected in assignments:
            mask = np.zeros(len(pooled), dtype=bool)
            mask[list(selected)] = True
            delta = float(pooled[mask].mean() - pooled[~mask].mean())
            extreme += abs(delta) >= abs(difference) - 1e-12
        result["exact_two_sided_permutation_p"] = extreme / len(assignments)
        result["permutation_assignments"] = len(assignments)
        result["permutation_limitation"] = "assumes exchangeability; condition/processing confounding is unresolved; no multiplicity correction"
    return result


def fit_panels(matrix, names, samples, input_sha256):
    """Only proliferating rows select controls, independently of held-out values."""
    train = [i for i, sample in enumerate(samples) if sample["condition"] == "proliferating"]
    training_ids = [samples[i]["accession"] for i in train]
    if len(train) < 2:
        raise ScoreError("at least two proliferating baseline samples required")
    fits, refusals = {}, {}
    for key, entry in GENE_SETS.items():
        try:
            fits[key] = fit_module_score(
                matrix, names, entry["symbols"], train, seed=0, n_ctrl=5, n_bins=20,
                block_from_controls=set(ORTHOGONAL), orthogonal=entry["orthogonal"],
                citation=entry["citation"], gene_set_key=key,
                gene_set_source_status=entry["source_status"], training_sample_ids=training_ids,
                training_input_sha256=input_sha256,
            )
        except ScoreError as exc:
            refusals[key] = str(exc)
    try:
        fits["fridman_signed"] = fit_signed_fridman_score(
            matrix, names, train, seed=0, n_ctrl=5, n_bins=20,
            training_sample_ids=training_ids, training_input_sha256=input_sha256,
        )
    except ScoreError as exc:
        refusals["fridman_signed"] = str(exc)
    return fits, refusals, train


def evaluate(raw_dir: Path, snapshot_dir: Path, output_dir: Path, download=False):
    manifest = json.loads((snapshot_dir / "sources.json").read_text(encoding="utf-8"))
    raw_dir.mkdir(parents=True, exist_ok=True)
    counts_path = raw_dir / "counts.txt.gz"
    if download and not counts_path.exists():
        with urllib.request.urlopen(manifest["counts"]["url"], timeout=90) as response:
            counts_path.write_bytes(response.read())
    if not counts_path.exists() or sha256(counts_path) != manifest["counts"]["sha256"]:
        raise ScoreError("actual GEO count input missing or hash mismatch; use --download or restore pinned input")
    for name in ("samples.json", "hgnc_mapping.tsv.gz", "PLAN.md"):
        if sha256(snapshot_dir / name) != manifest["snapshots"][name]:
            raise ScoreError(f"pinned source/plan snapshot hash mismatch: {name}")
    samples = json.loads((snapshot_dir / "samples.json").read_text(encoding="utf-8"))
    with gzip.open(snapshot_dir / "hgnc_mapping.tsv.gz", "rt", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    mapping = {row["ensembl_gene_id"]: row["symbol"] for row in rows}
    if len(mapping) != len(rows) or any(not ENSEMBL.fullmatch(identifier) for identifier in mapping):
        raise ScoreError("invalid/duplicate pinned mapping identifier")
    matrix, raw, names, samples, statistics = read_counts(counts_path, mapping, samples)
    groups = {condition: [i for i, sample in enumerate(samples) if sample["condition"] == condition] for condition in CONDITIONS}
    if any(len(rows) != 3 for rows in groups.values()):
        raise ScoreError("prespecified GSE268487 run requires exactly three libraries per condition")
    output_dir.mkdir(parents=True, exist_ok=True)
    fits, refusals, train = fit_panels(matrix, names, samples, manifest["counts"]["sha256"])
    scores, panels = {}, {}
    def summarize(vector):
        return {
            "scores_by_accession": {sample["accession"]: float(value) for sample, value in zip(samples, vector)},
            "means_by_condition": {condition: float(np.asarray(vector)[rows].mean()) for condition, rows in groups.items()},
            "senescent_vs_quiescent_held_out": descriptive_contrast(vector, groups["senescent"], groups["quiescent"], permutation=True),
            "senescent_vs_proliferating_calibration": descriptive_contrast(vector, groups["senescent"], train),
            "quiescent_vs_proliferating_calibration": descriptive_contrast(vector, groups["quiescent"], train),
        }
    for key, fit in fits.items():
        _write_json(output_dir / f"{key}_fit.json", fit)
        if key == "fridman_signed":
            transformed = transform_signed_fridman_score(matrix, names, fit)
            coverage = {"up": transformed["up_coverage"], "down": transformed["down_coverage"]}
            extra = {"up_missing": transformed["up_missing"], "down_missing": transformed["down_missing"],
                     "up_component": summarize(transformed["up_scores"]),
                     "down_component": summarize(transformed["down_scores"]),
                     "source_status": "published_gene_sets_project_signed_method"}
        else:
            transformed = transform_module_score(matrix, names, fit)
            coverage = transformed["coverage"]
            present = [names.index(gene) for gene in fit["signature"] if gene in names]
            extra = {"missing": transformed["missing"], "n_present": len(present),
                     "source_status": GENE_SETS[key]["source_status"],
                     "n_genes_expressed_by_condition": {condition: int((raw[rows][:, present].sum(axis=0) > 0).sum()) for condition, rows in groups.items()}}
        scores[key] = transformed["scores"]
        panels[key] = {"coverage": coverage, **extra, **summarize(transformed["scores"])}
    markers = {}
    for gene in ORTHOGONAL:
        if gene in names:
            values = matrix[:, names.index(gene)]
            reference = values[train]
            std = float(reference.std()) or 1.0
            markers[gene] = {"logCPM": summarize(values), "baseline_z": summarize((values - reference.mean()) / std),
                             "overlapping_panel_keys": [key for key, entry in GENE_SETS.items() if gene in entry["symbols"]]}
    forbidden = set(ORTHOGONAL).union(*(set(entry["symbols"]) for entry in GENE_SETS.values()))
    random_results = []
    for seed in range(1000, 1100):
        signature = random_signature(names, len(SENMAYO_HUMAN), seed, forbidden)
        fit = fit_module_score(matrix, names, signature, train, seed=0, n_ctrl=5, n_bins=20,
                               block_from_controls=forbidden, orthogonal=(), citation="Uniform random control, prespecified seeds 1000..1099",
                               gene_set_source_status="custom_unverified")
        vector = transform_module_score(matrix, names, fit)["scores"]
        random_results.append({"seed": seed, "signature": list(signature), **summarize(vector)})
    _write_json(output_dir / "random_controls.json", random_results)
    random_deltas = [entry["senescent_vs_quiescent_held_out"]["mean_difference"] for entry in random_results]
    random_aucs = [entry["senescent_vs_quiescent_held_out"]["directional_auroc"] for entry in random_results]
    empirical_rank = None
    if "senmayo" in panels:
        target = abs(panels["senmayo"]["senescent_vs_quiescent_held_out"]["mean_difference"])
        empirical_rank = (1 + sum(abs(delta) >= target for delta in random_deltas)) / (1 + len(random_deltas))
    overlaps = {f"{a}__{b}": sorted(set(GENE_SETS[a]["symbols"]) & set(GENE_SETS[b]["symbols"]))
                for a, b in itertools.combinations(GENE_SETS, 2)}
    correlations = {f"{a}__{b}": pearson_correlation(scores[a], scores[b]) for a, b in itertools.combinations(scores, 2)}
    summary = {
        "accession": ACCESSION, "sources": manifest, "software": {"numpy": np.__version__,
            "importer_version": IMPORTER_VERSION, "importer_source_sha256": sha256(Path(__file__)),
            "importer_correction_receipt_sha256": sha256(snapshot_dir / "IMPORTER_CORRECTION.md")},
        "import": statistics, "samples_in_matrix_order": samples,
        "training_sample_ids": [samples[i]["accession"] for i in train],
        "held_out_sample_ids": [sample["accession"] for i, sample in enumerate(samples) if i not in train],
        "split_basis": "condition baseline only; donor and experiment blocks unavailable; no donor independence claim",
        "panels": panels, "refusals": refusals, "orthogonal_expression_markers": markers,
        "signature_overlaps": overlaps, "cross_score_pearson_all_9_libraries": correlations,
        "random_controls": {"n": 100, "signature_size": 125, "seeds": [1000, 1099],
                            "mean_difference_quantiles": np.quantile(random_deltas, [0, .025, .5, .975, 1]).tolist(),
                            "directional_auroc_quantiles": np.quantile(random_aucs, [0, .025, .5, .975, 1]).tolist(),
                            "n_directional_auroc_one": sum(auc == 1 for auc in random_aucs),
                            "n_directional_auroc_zero": sum(auc == 0 for auc in random_aucs),
                            "senmayo_absolute_delta_empirical_rank": empirical_rank,
                            "limitation": "uniform, not expression matched; rank is not a validated specificity test"},
        "functional_assay_linkage": "not available; separate IMR90 qPCR samples were not joined",
        "interpretation_limits": ["single LF1 cell line", "donor metadata unavailable", "reported replicates do not establish independent donors",
                                  "condition confounded with sample processing and lanes", "linked publication is a preprint",
                                  "no tissue generalization or rejuvenation result", "no patient age or diagnostic inference"],
    }
    _write_json(output_dir / "summary.json", summary)
    with (output_dir / "expression_logCPM.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample_id", *names])
        for sample, row in zip(samples, matrix):
            writer.writerow([sample["accession"], *row])
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, default=Path("data/public/GSE268487"))
    parser.add_argument("--snapshot-dir", type=Path, default=Path("validation/GSE268487"))
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/GSE268487"))
    parser.add_argument("--download", action="store_true")
    args = parser.parse_args(argv)
    try:
        summary = evaluate(args.raw_dir, args.snapshot_dir, args.output_dir, args.download)
    except (ScoreError, OSError, ValueError) as exc:
        parser.error(str(exc))
    print(json.dumps({"accession": summary["accession"], "summary": str(args.output_dir / "summary.json"),
                      "refusals": summary["refusals"]}, indent=2))


if __name__ == "__main__":
    main()

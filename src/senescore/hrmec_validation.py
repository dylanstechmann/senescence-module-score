"""Prespecified external check of the control-bin senescence scores on GEO GSE160356.

GSE160356 holds nine RNA-seq libraries from human retinal microvascular endothelial cells
(HRMEC): three early-passage, three late-passage (replicative) and three etoposide-treated.
The question and every choice below are fixed in ``validation/GSE160356/PLAN.md`` before any
score is calculated; ``evaluate`` refuses to run if that file no longer matches the hash pinned
in ``sources.json``.

What is different from the earlier GSE268487 check, and why:

* The held-out libraries (late passage, etoposide) never influence control selection.
* The early-passage libraries are the baseline, so scoring them with a fit that used them would
  centre their scores near zero by construction. The primary protocol therefore scores each
  early-passage library with a fit made from the *other two* (leave-one-library-out). The all-three
  fit is kept as a labelled sensitivity analysis.
* A plain-language verdict rule (pass / ranking only / fail) is fixed in advance, together with a
  random-panel comparison under the same protocol.

This importer never infers clones, donors or pairing from replicate numbers, never changes a
signature, and never picks seeds after seeing results.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import io
import itertools
import json
import re
import subprocess
import tarfile
import urllib.request
from pathlib import Path

import numpy as np

from senescore.genes import GENE_SETS, ORTHOGONAL, SENMAYO_HUMAN
from senescore.public_validation import descriptive_contrast, read_counts
from senescore.score import (
    ScoreError, fit_module_score, fit_signed_fridman_score, pearson_correlation, random_signature,
    transform_module_score, transform_signed_fridman_score,
)

ACCESSION = "GSE160356"
EVALUATOR_VERSION = 1
BASELINE = "early_passage"
INDUCED = ("late_passage", "etoposide")
CONDITIONS = (BASELINE, *INDUCED)
TITLE = re.compile(r"^(EP_CT|LP_CT|ETO)_([1-3])$")
PREFIX = {"EP_CT": ("early_passage", "ep"), "LP_CT": ("late_passage", "lp"), "ETO": ("etoposide", "eto")}
TREATMENT = {"early_passage": "Early Passage", "late_passage": "Late Passage", "etoposide": "Etoposide Treatment"}
MEMBER = re.compile(r"^(GSM\d+)_[A-Za-z0-9_.\-]+\.genename\.htcounts\.txt\.gz$")
MAX_MEMBER_BYTES = 5_000_000
RANDOM_SEEDS = range(1000, 1100)
RANK_THRESHOLD = 0.05
SIGNATURE_SIZE = len(SENMAYO_HUMAN)
CONTRASTS = (("late_passage_vs_early_passage", "late_passage", BASELINE),
             ("etoposide_vs_early_passage", "etoposide", BASELINE),
             ("etoposide_vs_late_passage", "etoposide", "late_passage"))
PRIMARY_CONTRASTS = ("late_passage_vs_early_passage", "etoposide_vs_early_passage")
VERDICT_RULE = (
    "Applies to SenMayo only, under the primary protocol. 'pass': directional AUROC is 1.0 for both "
    "late-passage-vs-early-passage and etoposide-vs-early-passage AND SenMayo's absolute mean difference "
    f"has a two-sided empirical rank <= {RANK_THRESHOLD} against {len(RANDOM_SEEDS)} uniform random "
    f"{SIGNATURE_SIZE}-gene panels for both contrasts. 'ranking_only_not_distinguishable_from_random_panels': "
    "both AUROCs are 1.0 but at least one rank exceeds the threshold. 'fail': any AUROC is below 1.0. "
    "'not_scored': SenMayo was refused for low gene coverage."
)
PRIMARY_METRIC = "directional_auroc_vs_early_passage"
RANDOM_BASELINE = "random_uniform_panels_100"


def sha256_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(Path(path).read_bytes())


def write_json(path: Path, value) -> None:
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")


# --------------------------------------------------------------------------- metadata

def parse_geo_samples(text: str) -> dict:
    """Extract exact accession, title and condition facts; absent clones and donors stay absent."""
    series = {}
    for key in ("!Series_title", "!Series_overall_design", "!Series_web_link", "!Series_status",
                "!Series_submission_date"):
        match = re.search(rf"^{re.escape(key)} = (.*)$", text, flags=re.M)
        if match:
            series[key.removeprefix("!Series_")] = match.group(1)
    blocks = re.split(r"^\^SAMPLE = ", text, flags=re.M)[1:]
    samples = []
    for block in blocks:
        accession, _, body = block.partition("\n")
        fields: dict[str, list[str]] = {}
        for line in body.splitlines():
            if line.startswith("^"):
                break
            if line.startswith("!Sample_") and " = " in line:
                key, value = line.split(" = ", 1)
                fields.setdefault(key, []).append(value)
        titles = fields.get("!Sample_title", [])
        match = TITLE.fullmatch(titles[0]) if len(titles) == 1 else None
        if not match or fields.get("!Sample_organism_ch1") != ["Homo sapiens"]:
            raise ScoreError("unexpected or incomplete GSE160356 sample metadata")
        condition, short = PREFIX[match.group(1)]
        characteristics = fields.get("!Sample_characteristics_ch1", [])
        if (f"treatment: {TREATMENT[condition]}" not in characteristics
                or "tissue: Endothelial Cells" not in characteristics
                or fields.get("!Sample_source_name_ch1") != ["Human Retinal Microvascular Endothelial Cells"]):
            raise ScoreError("sample metadata has an unexpected cell type or treatment")
        files = fields.get("!Sample_supplementary_file_1", [])
        member = files[0].rsplit("/", 1)[-1] if len(files) == 1 else ""
        member_match = MEMBER.fullmatch(member)
        if not member_match or member_match.group(1) != accession.strip():
            raise ScoreError("sample supplementary count file does not belong to its accession")
        passage = next((item.removeprefix("passage: ") for item in characteristics if item.startswith("passage: ")), None)
        samples.append({
            "accession": accession.strip(), "title": titles[0], "condition": condition,
            "reported_replicate": int(match.group(2)), "column": f"{short}{match.group(2)}",
            "cell_type": "HRMEC", "reported_passage": passage, "count_file": member,
            "donor_id": None, "clone_id": None, "experiment_block": None,
            "source_fields": {key: values for key, values in fields.items() if key in (
                "!Sample_organism_ch1", "!Sample_characteristics_ch1", "!Sample_source_name_ch1",
                "!Sample_treatment_protocol_ch1", "!Sample_extract_protocol_ch1",
                "!Sample_data_processing", "!Sample_status", "!Sample_platform_id",
                "!Sample_instrument_model", "!Sample_library_strategy")},
        })
    if not samples:
        raise ScoreError("GEO metadata contains no samples")
    if (len({s["column"] for s in samples}) != len(samples) or len({s["accession"] for s in samples}) != len(samples)
            or len({s["count_file"] for s in samples}) != len(samples)):
        raise ScoreError("duplicate sample accession, condition/replicate or count file")
    return {"series": series, "samples": sorted(samples, key=lambda s: s["accession"])}


# --------------------------------------------------------------------------- inputs

def extract_members(tar_path: Path, expected: dict[str, dict]) -> dict[str, bytes]:
    """Read exactly the pinned members of the GEO RAW archive into memory, nothing else."""
    members: dict[str, bytes] = {}
    with tarfile.open(tar_path, "r:") as archive:
        for info in archive.getmembers():
            if info.name not in expected:
                raise ScoreError(f"unexpected member in the GEO archive: {info.name!r}")
            if not info.isfile() or info.size > MAX_MEMBER_BYTES or info.name in members:
                raise ScoreError(f"archive member is not a single bounded regular file: {info.name!r}")
            handle = archive.extractfile(info)
            raw = handle.read(MAX_MEMBER_BYTES + 1) if handle else b""
            if len(raw) != info.size or sha256_bytes(raw) != expected[info.name]["sha256"]:
                raise ScoreError(f"archive member does not match its pinned hash: {info.name}")
            members[info.name] = raw
    if set(members) != set(expected):
        raise ScoreError("the GEO archive is missing pinned members")
    return members


def assemble_counts(members: dict[str, bytes], samples: list[dict]) -> bytes:
    """Join per-library HTSeq tables into one deterministic gene-by-library TSV (gzip, mtime 0)."""
    ordered = sorted(samples, key=lambda s: s["column"])
    tables = []
    for sample in ordered:
        raw = members[sample["count_file"]]
        try:
            lines = gzip.decompress(raw).decode("utf-8").splitlines()
        except (OSError, EOFError, UnicodeDecodeError) as exc:
            raise ScoreError(f"unreadable count file for {sample['accession']}") from exc
        rows = [line.split("\t") for line in lines]
        if not rows or any(len(row) != 2 or not row[0].strip() for row in rows):
            raise ScoreError(f"count file for {sample['accession']} is not a two-column table")
        tables.append(rows)
    identifiers = [row[0] for row in tables[0]]
    if len(identifiers) != len(set(identifiers)):
        raise ScoreError("duplicate identifier in a count file")
    for sample, rows in zip(ordered, tables):
        if [row[0] for row in rows] != identifiers:
            raise ScoreError(f"identifiers of {sample['accession']} differ from the other libraries")
    buffer = io.BytesIO()
    with gzip.GzipFile(fileobj=buffer, mode="wb", mtime=0) as handle:
        handle.write(("\t".join(["gene", *[s["column"] for s in ordered]]) + "\n").encode("utf-8"))
        for position, identifier in enumerate(identifiers):
            handle.write(("\t".join([identifier, *[table[position][1] for table in tables]]) + "\n").encode("utf-8"))
    return buffer.getvalue()


def git_state(repo_root: Path) -> dict:
    """The checked-out revision and whether tracked files differ from it, or unknowns."""
    try:
        head = subprocess.run(["git", "-C", str(repo_root), "rev-parse", "HEAD"], capture_output=True,
                              text=True, timeout=15, check=False)
        status = subprocess.run(["git", "-C", str(repo_root), "status", "--porcelain", "--untracked-files=no"],
                                capture_output=True, text=True, timeout=15, check=False)
    except (OSError, subprocess.SubprocessError):
        return {"revision": None, "tracked_changes": None}
    revision = head.stdout.strip() if head.returncode == 0 and re.fullmatch(r"[0-9a-f]{40}", head.stdout.strip()) else None
    changes = bool(status.stdout.strip()) if status.returncode == 0 else None
    return {"revision": revision, "tracked_changes": changes}


# --------------------------------------------------------------------------- scoring

def _fit_for(key, matrix, names, samples, rows, input_sha256):
    ids = [samples[i]["accession"] for i in rows]
    if key == "fridman_signed":
        return fit_signed_fridman_score(matrix, names, rows, seed=0, n_ctrl=5, n_bins=20,
                                        training_sample_ids=ids, training_input_sha256=input_sha256)
    entry = GENE_SETS[key]
    return fit_module_score(matrix, names, entry["symbols"], rows, seed=0, n_ctrl=5, n_bins=20,
                            block_from_controls=set(ORTHOGONAL), orthogonal=entry["orthogonal"],
                            citation=entry["citation"], gene_set_key=key,
                            gene_set_source_status=entry["source_status"], training_sample_ids=ids,
                            training_input_sha256=input_sha256)


def _transform(key, matrix, names, fit):
    if key == "fridman_signed":
        return np.asarray(transform_signed_fridman_score(matrix, names, fit)["scores"], dtype=float)
    return np.asarray(transform_module_score(matrix, names, fit)["scores"], dtype=float)


def two_protocol_scores(baseline_rows, fit_fn, transform_fn):
    """Return ``(primary, sensitivity, full_fit)``: score vectors over every library, plus the all-baseline fit.

    Sensitivity: every library is scored with the fit from all baseline libraries. Primary: each
    baseline library is scored with the fit from the *other* baseline libraries, so no baseline
    library helps choose the controls used to score it; non-baseline libraries still use the
    all-baseline fit.
    """
    full_fit = fit_fn(list(baseline_rows))
    sensitivity = transform_fn(full_fit)
    primary = np.array(sensitivity, dtype=float, copy=True)
    for held in baseline_rows:
        train = [row for row in baseline_rows if row != held]
        primary[held] = transform_fn(fit_fn(train))[held]
    return primary, np.asarray(sensitivity, dtype=float), full_fit


def contrasts_for(vector, groups, *, permutation: bool):
    return {name: descriptive_contrast(vector, groups[positive], groups[negative], permutation=permutation)
            for name, positive, negative in CONTRASTS}


def empirical_rank(target: float, null_values) -> float:
    """Two-sided empirical rank of ``abs(target)`` among ``abs`` null values (includes the target)."""
    return (1 + sum(abs(value) >= abs(target) - 1e-12 for value in null_values)) / (1 + len(null_values))


def verdict(senmayo_primary, ranks, refused: bool) -> str:
    if refused or senmayo_primary is None:
        return "not_scored"
    if all(senmayo_primary[name]["directional_auroc"] == 1.0 for name in PRIMARY_CONTRASTS):
        if all(ranks[name] is not None and ranks[name] <= RANK_THRESHOLD for name in PRIMARY_CONTRASTS):
            return "pass"
        return "ranking_only_not_distinguishable_from_random_panels"
    return "fail"


def build_receipt(*, created_utc, state, input_hashes, baseline_ids, held_out_ids, panel_keys, metrics, leakage):
    return {
        "schema": "regen-workbench/evaluation-receipt/1",
        "producer": {"tool": "senescore.hrmec_validation", "version": str(EVALUATOR_VERSION),
                     "code_revision": state["revision"]},
        "created_utc": created_utc,
        "input_sha256": sorted(set(input_hashes)),
        "folds": [{"train_groups": sorted(baseline_ids), "test_groups": sorted(held_out_ids)}],
        "model_names": sorted([*panel_keys, RANDOM_BASELINE]),
        "metric_names": sorted(metrics),
        "reported_leakage": leakage,
        "seed": 0,
        "tree_dirty": state["tracked_changes"],
        "baseline_scoring": {
            "protocol": "each early-passage library is scored with a fit from the other two (leave-one-out); "
                        "held-out libraries use the fit from all three early-passage libraries",
            "note": "These leave-one-out baseline scores are development-side computations; no sealed library is a "
                    "test group in them, so they are not listed as folds."},
        "disclaimer": "Procedural receipt for a descriptive 3-vs-3 check. It does not establish biological "
                      "senescence, independence of the libraries, or anything about tissue.",
    }


# --------------------------------------------------------------------------- run

def load_inputs(raw_dir: Path, snapshot_dir: Path, hgnc_path: Path, *, download: bool):
    manifest = json.loads((snapshot_dir / "sources.json").read_text(encoding="utf-8"))
    raw_dir.mkdir(parents=True, exist_ok=True)
    tar_path = raw_dir / "GSE160356_RAW.tar"
    if download and not tar_path.exists():
        with urllib.request.urlopen(manifest["raw_archive"]["url"], timeout=120) as response:
            tar_path.write_bytes(response.read(manifest["raw_archive"]["bytes"] + 1))
    if not tar_path.exists() or sha256_file(tar_path) != manifest["raw_archive"]["sha256"]:
        raise ScoreError("actual GEO count archive missing or hash mismatch; use --download or restore pinned input")
    for name in ("samples.json", "PLAN.md"):
        if sha256_file(snapshot_dir / name) != manifest["snapshots"][name]:
            raise ScoreError(f"pinned source/plan snapshot hash mismatch: {name}")
    if sha256_file(hgnc_path) != manifest["hgnc"]["sha256"]:
        raise ScoreError("pinned HGNC mapping snapshot hash mismatch")
    snapshot = json.loads((snapshot_dir / "samples.json").read_text(encoding="utf-8"))
    samples = snapshot["samples"]
    members = extract_members(tar_path, manifest["members"])
    with gzip.open(hgnc_path, "rt", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    mapping = {row["ensembl_gene_id"]: row["symbol"] for row in rows}
    if len(mapping) != len(rows):
        raise ScoreError("duplicate identifier in the pinned HGNC mapping")
    return manifest, samples, members, mapping


def evaluate(raw_dir: Path, snapshot_dir: Path, output_dir: Path, hgnc_path: Path, *, download=False,
             validate_only=False, allow_dirty=False, repo_root=None, now=None):
    manifest, samples, members, mapping = load_inputs(raw_dir, snapshot_dir, hgnc_path, download=download)
    counts_bytes = assemble_counts(members, samples)
    counts_sha256 = sha256_bytes(counts_bytes)
    counts_path = raw_dir / "counts_combined.tsv.gz"
    counts_path.write_bytes(counts_bytes)
    matrix, raw, names, ordered, statistics = read_counts(counts_path, mapping, samples)
    groups = {c: [i for i, s in enumerate(ordered) if s["condition"] == c] for c in CONDITIONS}
    if any(len(rows) != 3 for rows in groups.values()):
        raise ScoreError("prespecified GSE160356 run requires exactly three libraries per condition")
    if validate_only:
        return {"accession": ACCESSION, "validated": True, "libraries": len(ordered),
                "ensembl_rows": statistics["ensembl_rows"], "mapped_rows": statistics["mapped_rows"],
                "mapped_symbols": statistics["mapped_symbols"], "combined_counts_sha256": counts_sha256}

    state = git_state(Path(repo_root) if repo_root else Path(__file__).resolve().parents[2])
    if state["tracked_changes"] and not allow_dirty:
        raise ScoreError("tracked files differ from the committed revision; commit or stash before the "
                         "prespecified run so the recorded code revision describes what ran")
    created_utc = (now or dt.datetime.now(dt.timezone.utc)).isoformat(timespec="seconds")
    output_dir.mkdir(parents=True, exist_ok=True)
    baseline = groups[BASELINE]
    input_hashes = [info["sha256"] for info in manifest["members"].values()] + [manifest["hgnc"]["sha256"]]

    trained_rows: set[int] = set()  # every library any fit was allowed to see

    def logged(factory):
        def run(rows):
            trained_rows.update(rows)
            return factory(rows)
        return run

    panel_keys = [*GENE_SETS, "fridman_signed"]
    panels, refusals, score_primary, full_fits = {}, {}, {}, {}
    for key in panel_keys:
        try:
            primary, sensitivity, full_fit = two_protocol_scores(
                baseline,
                logged(lambda rows, key=key: _fit_for(key, matrix, names, ordered, rows, counts_sha256)),
                lambda fit, key=key: _transform(key, matrix, names, fit))
        except ScoreError as exc:
            refusals[key] = str(exc)
            continue
        full_fits[key] = full_fit
        score_primary[key] = primary
        coverage = (full_fit["up_fit"]["coverage"], full_fit["down_fit"]["coverage"]) if key == "fridman_signed" \
            else full_fit["coverage"]
        panels[key] = {
            "coverage": coverage,
            "source_status": "published_gene_sets_project_signed_method" if key == "fridman_signed"
            else GENE_SETS[key]["source_status"],
            "primary": {"scores_by_accession": {s["accession"]: float(v) for s, v in zip(ordered, primary)},
                        "contrasts": contrasts_for(primary, groups, permutation=True)},
            "sensitivity_all_baseline_fit": {
                "scores_by_accession": {s["accession"]: float(v) for s, v in zip(ordered, sensitivity)},
                "contrasts": contrasts_for(sensitivity, groups, permutation=True)},
        }

    forbidden = set(ORTHOGONAL).union(*(set(entry["symbols"]) for entry in GENE_SETS.values()))
    random_rows, forbidden_in_random_panels = [], 0
    for seed in RANDOM_SEEDS:
        signature = random_signature(names, SIGNATURE_SIZE, seed, forbidden)
        forbidden_in_random_panels += len(set(signature) & forbidden)
        primary, _sensitivity, _fit = two_protocol_scores(
            baseline,
            logged(lambda rows, signature=signature: fit_module_score(
                matrix, names, signature, rows, seed=0, n_ctrl=5, n_bins=20, block_from_controls=forbidden,
                orthogonal=(), citation="Uniform random control, prespecified seeds 1000..1099",
                gene_set_source_status="custom_unverified")),
            lambda fit: np.asarray(transform_module_score(matrix, names, fit)["scores"], dtype=float))
        row = {"seed": seed}
        for name, positive, negative in CONTRASTS:
            result = descriptive_contrast(primary, groups[positive], groups[negative])
            row[f"{name}__mean_difference"] = result["mean_difference"]
            row[f"{name}__directional_auroc"] = result["directional_auroc"]
        random_rows.append(row)
    with (snapshot_dir / "random_scores.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(random_rows[0]))
        writer.writeheader()
        writer.writerows(random_rows)

    senmayo = panels.get("senmayo")
    ranks = {}
    random_summary = {"n": len(random_rows), "signature_size": SIGNATURE_SIZE,
                      "seeds": [RANDOM_SEEDS.start, RANDOM_SEEDS.stop - 1],
                      "excluded_from_signatures_and_controls": "union of all registered panels plus CDKN1A/CDKN2A",
                      "limitation": "uniform, not expression matched; the rank is not a validated specificity test",
                      "contrasts": {}}
    for name, _positive, _negative in CONTRASTS:
        deltas = [row[f"{name}__mean_difference"] for row in random_rows]
        aucs = [row[f"{name}__directional_auroc"] for row in random_rows]
        target = senmayo["primary"]["contrasts"][name]["mean_difference"] if senmayo else None
        ranks[name] = empirical_rank(target, deltas) if senmayo else None
        random_summary["contrasts"][name] = {
            "mean_difference_quantiles_0_2.5_50_97.5_100": np.quantile(deltas, [0, .025, .5, .975, 1]).tolist(),
            "directional_auroc_quantiles_0_2.5_50_97.5_100": np.quantile(aucs, [0, .025, .5, .975, 1]).tolist(),
            "n_directional_auroc_one": int(sum(a == 1.0 for a in aucs)),
            "n_directional_auroc_zero": int(sum(a == 0.0 for a in aucs)),
            "senmayo_absolute_delta_empirical_rank": ranks[name]}

    result_verdict = verdict(senmayo["primary"]["contrasts"] if senmayo else None, ranks, "senmayo" in refusals)
    markers = {}
    for gene in ORTHOGONAL:
        if gene in names:
            values = matrix[:, names.index(gene)]
            reference = values[baseline]
            std = float(reference.std()) or 1.0
            markers[gene] = {"logCPM_by_accession": {s["accession"]: float(v) for s, v in zip(ordered, values)},
                             "baseline_z_all_baseline_fit": {s["accession"]: float((v - reference.mean()) / std)
                                                             for s, v in zip(ordered, values)},
                             "overlapping_panel_keys": [k for k, e in GENE_SETS.items() if gene in e["symbols"]]}
    overlaps = {f"{a}__{b}": sorted(set(GENE_SETS[a]["symbols"]) & set(GENE_SETS[b]["symbols"]))
                for a, b in itertools.combinations(GENE_SETS, 2)}
    correlations = {f"{a}__{b}": pearson_correlation(score_primary[a], score_primary[b])
                    for a, b in itertools.combinations(score_primary, 2)}
    baseline_ids = [ordered[i]["accession"] for i in baseline]
    held_out_ids = [ordered[i]["accession"] for i in range(len(ordered)) if i not in baseline]
    leakage = [{"id": "held_out_libraries_used_to_fit_controls", "count": len(trained_rows - set(baseline))},
               {"id": "registered_panel_genes_in_random_panels", "count": forbidden_in_random_panels}]
    results = {
        "accession": ACCESSION, "evaluator_version": EVALUATOR_VERSION, "created_utc": created_utc,
        "code_revision": state["revision"], "tracked_tree_changes": state["tracked_changes"],
        "plan_sha256": manifest["snapshots"]["PLAN.md"], "combined_counts_sha256": counts_sha256,
        "software": {"numpy": np.__version__},
        "import": statistics, "libraries_in_matrix_order": [s["accession"] for s in ordered],
        "split": {"baseline_libraries": baseline_ids, "held_out_libraries": held_out_ids,
                  "leave_one_out_baseline_folds": [
                      {"train": [ordered[j]["accession"] for j in baseline if j != i], "scored": ordered[i]["accession"]}
                      for i in baseline],
                  "basis": "condition baseline only; clone, donor and experiment block are not deposited; no independence claim"},
        "panels": panels, "refusals": refusals, "orthogonal_expression_markers": markers,
        "signature_overlaps": overlaps, "cross_score_pearson_primary_all_9_libraries": correlations,
        "random_controls": random_summary,
        "verdict": {"rule": VERDICT_RULE, "value": result_verdict, "senmayo_empirical_ranks": ranks,
                    "senmayo_primary_directional_auroc": {
                        name: senmayo["primary"]["contrasts"][name]["directional_auroc"] for name in PRIMARY_CONTRASTS
                    } if senmayo else None},
        "interpretation_limits": [
            "three libraries per condition; the smallest possible exact two-sided permutation p is 0.1",
            "series text says three independent clones, but no sample-level clone or pairing is deposited; none is inferred",
            "passage ranges differ by condition, so passage is confounded with condition",
            "no per-library functional senescence assay is deposited with the series",
            "the linked study derived its own endothelial signature (EndoSEN) from related data; it is not scored here",
            "uniform random panels are not expression matched and do not establish specificity",
            "no tissue generalization, rejuvenation or diagnostic claim"],
    }
    write_json(snapshot_dir / "results.json", results)
    metrics = [PRIMARY_METRIC, "mean_difference_vs_early_passage", "cohen_d_vs_early_passage",
               "exact_two_sided_permutation_p"]
    receipt = build_receipt(created_utc=created_utc, state=state, input_hashes=input_hashes,
                            baseline_ids=baseline_ids, held_out_ids=held_out_ids,
                            panel_keys=list(panels), metrics=metrics, leakage=leakage)
    write_json(snapshot_dir / "evaluation_receipt.json", receipt)
    for key, fit in full_fits.items():
        write_json(output_dir / f"{key}_fit.json", fit)
    with (output_dir / "expression_logCPM.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample_id", *names])
        for sample, row in zip(ordered, matrix):
            writer.writerow([sample["accession"], *row])
    return results


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--raw-dir", type=Path, default=Path("data/public/GSE160356"))
    parser.add_argument("--snapshot-dir", type=Path, default=Path("validation/GSE160356"))
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/GSE160356"))
    parser.add_argument("--hgnc-mapping", type=Path, default=Path("validation/GSE268487/hgnc_mapping.tsv.gz"))
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--validate-only", action="store_true",
                        help="check hashes and table structure and print counts only; calculates no score")
    parser.add_argument("--allow-dirty-tree", action="store_true",
                        help="development only: a prespecified run should come from a clean committed tree")
    args = parser.parse_args(argv)
    try:
        summary = evaluate(args.raw_dir, args.snapshot_dir, args.output_dir, args.hgnc_mapping,
                           download=args.download, validate_only=args.validate_only,
                           allow_dirty=args.allow_dirty_tree)
    except (ScoreError, OSError, ValueError, tarfile.TarError) as exc:
        parser.error(str(exc))
    if args.validate_only:
        print(json.dumps(summary, indent=2))
    else:
        print(json.dumps({"accession": summary["accession"], "verdict": summary["verdict"]["value"],
                          "results": str(args.snapshot_dir / "results.json"),
                          "refusals": summary["refusals"]}, indent=2))


if __name__ == "__main__":
    main()

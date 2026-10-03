"""Seurat-style control-gene module score, in numpy.

For each signature gene, subtract the average of control genes from the same
expression bin. A random gene set of the same size is the baseline, not a
second senescence clock.
"""

from __future__ import annotations

import hashlib
import json

import numpy as np

from senescore.genes import CITATION, ORTHOGONAL, SENMAYO_HUMAN


class ScoreError(ValueError):
    pass


FIT_ARTIFACT_FORMAT = "senescore-module-score-fit"
FIT_ARTIFACT_VERSION = 1


def _json_sha256(value) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def cohen_d(group_a: np.ndarray, group_b: np.ndarray) -> float:
    group_a, group_b = np.asarray(group_a, dtype=float), np.asarray(group_b, dtype=float)
    if any(g.ndim != 1 or not np.isfinite(g).all() for g in (group_a, group_b)):
        raise ScoreError("effect size requires finite one-dimensional groups")
    n1, n2 = len(group_a), len(group_b)
    if n1 < 2 or n2 < 2:
        raise ScoreError("need at least two samples per group")
    v1, v2 = group_a.var(ddof=1), group_b.var(ddof=1)
    pooled = np.sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
    if pooled == 0:
        if group_a.mean() == group_b.mean():
            return 0.0
        raise ScoreError("effect size is undefined: unequal constant groups have zero pooled variance")
    return float((group_a.mean() - group_b.mean()) / pooled)


def pearson_correlation(x: np.ndarray, y: np.ndarray) -> float | None:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or len(x) < 2 or len(x) != len(y):
        return None
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        return None
    vx = x - x.mean()
    vy = y - y.mean()
    denom = float(np.sqrt(np.sum(vx ** 2) * np.sum(vy ** 2)))
    if denom == 0.0:
        return None
    return float(np.sum(vx * vy) / denom)


def _average_ranks(values: np.ndarray) -> np.ndarray:
    """Return one-based average ranks, assigning tied values their midrank."""
    order = np.argsort(values, kind="mergesort")
    sorted_values = values[order]
    starts = np.r_[0, np.flatnonzero(sorted_values[1:] != sorted_values[:-1]) + 1]
    stops = np.r_[starts[1:], len(values)]
    ranks = np.empty(len(values), dtype=float)
    for start, stop in zip(starts, stops):
        ranks[order[start:stop]] = (start + 1 + stop) / 2.0
    return ranks


def spearman_correlation(x: np.ndarray, y: np.ndarray) -> float | None:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or len(x) < 2 or len(x) != len(y):
        return None
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        return None
    rx = _average_ranks(x)
    ry = _average_ranks(y)
    return pearson_correlation(rx, ry)



def random_signature(gene_names: list[str], size: int, seed: int, forbidden: set[str]) -> tuple[str, ...]:
    rng = np.random.default_rng(seed)
    pool = [gene for gene in gene_names if gene.upper() not in forbidden]
    pick = rng.choice(pool, size=size, replace=False)
    return tuple(str(gene) for gene in pick)


def _validated_names(matrix: np.ndarray, gene_names: list[str], signature, n_ctrl: int, n_bins: int):
    matrix = np.asarray(matrix, dtype=np.float64)
    names = [gene.strip().upper() for gene in gene_names]
    signature = [gene.strip().upper() for gene in signature]
    if matrix.ndim != 2 or not matrix.shape[0] or matrix.shape[1] != len(names):
        raise ScoreError("matrix must be samples by genes")
    if not np.isfinite(matrix).all():
        raise ScoreError("expression values must be finite")
    for label, symbols in [("gene names", names), ("signature", signature)]:
        if not symbols or any(not s for s in symbols) or len(symbols) != len(set(symbols)):
            raise ScoreError(f"{label} must be nonempty, unique gene symbols")
    if any(isinstance(n, bool) or not isinstance(n, int) or n <= 0 for n in (n_ctrl, n_bins)):
        raise ScoreError("n_ctrl and n_bins must be positive integers")
    return matrix, names, signature


def _select_controls(means, names, signature, *, seed, n_ctrl, n_bins, block_from_controls):
    index = {gene: i for i, gene in enumerate(names)}
    present = [gene for gene in signature if gene in index]
    missing = [gene for gene in signature if gene not in index]
    coverage = len(present) / len(signature)
    if coverage < 0.6:
        raise ScoreError(f"only {coverage:.0%} of the signature is in the table")
    order = np.lexsort((np.array(names), means))
    bins = np.empty(len(names), dtype=np.int64)
    cuts = np.array_split(order, min(n_bins, len(names)))
    for b, cols in enumerate(cuts):
        bins[cols] = b
    rng = np.random.default_rng(seed)
    blocked = {gene.upper() for gene in signature}
    if block_from_controls:
        blocked |= {gene.upper() for gene in block_from_controls}
    sig_cols = np.array([index[gene] for gene in present])
    blocked_cols = np.array([index[gene] for gene in names if gene in blocked], dtype=np.int64)
    control_cols = []
    for col in sig_cols:
        pool = np.where(bins == bins[col])[0]
        pool = pool[~np.isin(pool, blocked_cols)]
        if len(pool) < n_ctrl:
            neighbor = np.where(np.abs(bins - bins[col]) <= 1)[0]
            pool = neighbor[~np.isin(neighbor, blocked_cols)]
        if len(pool) == 0:
            raise ScoreError("no control genes available")
        pool = np.array(sorted(pool, key=lambda i: names[i]))
        take = min(n_ctrl, len(pool))
        control_cols.append(rng.choice(pool, size=take, replace=False))
    return present, missing, coverage, sig_cols, control_cols


def _apply_controls(
    matrix,
    names,
    present,
    sig_cols,
    control_cols,
    *,
    z_reference_rows=None,
    orthogonal_genes=ORTHOGONAL,
    orthogonal_reference=None,
):
    sig_mean = matrix[:, sig_cols].mean(axis=1)
    ctrl_mean = np.mean([matrix[:, cols].mean(axis=1) for cols in control_cols], axis=0)
    scores = sig_mean - ctrl_mean
    index = {gene: i for i, gene in enumerate(names)}
    orthogonal = {}
    for gene in orthogonal_genes:
        if gene in index:
            col = matrix[:, index[gene]]
            if orthogonal_reference is not None:
                mean, std = orthogonal_reference[gene]
            else:
                reference = col if z_reference_rows is None else col[z_reference_rows]
                mean, std = float(reference.mean()), float(reference.std()) or 1.0
            orthogonal[gene] = ((col - mean) / std).tolist()
    return scores, orthogonal


def _validate_train_rows(train_rows, n_rows):
    rows = np.asarray(train_rows)
    if rows.ndim != 1 or rows.dtype.kind not in "iu" or len(rows) < 2 or len(set(rows.tolist())) != len(rows):
        raise ScoreError("train_rows must be at least two unique sample indexes")
    if int(rows.min()) < 0 or int(rows.max()) >= n_rows:
        raise ScoreError("train_rows out of range")
    return rows


def fit_module_score(
    matrix: np.ndarray,
    gene_names: list[str],
    signature: tuple[str, ...] | list[str],
    train_rows,
    seed: int = 0,
    n_ctrl: int = 5,
    n_bins: int = 20,
    block_from_controls: set[str] | None = None,
    orthogonal: tuple[str, ...] = ORTHOGONAL,
    citation: str = CITATION,
    *,
    gene_set_key: str | None = None,
    training_sample_ids: list[str] | None = None,
    training_input_sha256: str | None = None,
):
    """Fit and export the frozen controls and marker references for later batches."""
    matrix, names, signature = _validated_names(matrix, gene_names, signature, n_ctrl, n_bins)
    rows = _validate_train_rows(train_rows, len(matrix))
    if training_sample_ids is not None:
        if len(training_sample_ids) != len(rows) or any(not sample for sample in training_sample_ids):
            raise ScoreError("training_sample_ids must match train_rows and be nonempty")
        if len(training_sample_ids) != len(set(training_sample_ids)):
            raise ScoreError("training_sample_ids must be unique")
    if training_input_sha256 is not None and (
        len(training_input_sha256) != 64 or any(c not in "0123456789abcdef" for c in training_input_sha256.lower())
    ):
        raise ScoreError("training_input_sha256 must be a SHA-256 hex digest")

    present, missing, coverage, _sig_cols, control_cols = _select_controls(
        matrix[rows].mean(axis=0), names, signature, seed=seed, n_ctrl=n_ctrl, n_bins=n_bins,
        block_from_controls=block_from_controls,
    )
    index = {gene: i for i, gene in enumerate(names)}
    control_map = {
        gene: [names[i] for i in cols]
        for gene, cols in zip(present, control_cols)
    }
    orthogonal_reference = {}
    for gene in orthogonal:
        if gene in index:
            values = matrix[rows, index[gene]]
            std = float(values.std()) or 1.0
            orthogonal_reference[gene] = [float(values.mean()), std]

    artifact = {
        "format": FIT_ARTIFACT_FORMAT,
        "version": FIT_ARTIFACT_VERSION,
        "gene_set_key": gene_set_key,
        "signature": list(signature),
        "signature_sha256": _json_sha256(list(signature)),
        "feature_genes": list(names),
        "feature_schema_sha256": _json_sha256(list(names)),
        "training_rows": [int(i) for i in rows],
        "training_sample_ids": list(training_sample_ids) if training_sample_ids is not None else None,
        "training_sample_ids_sha256": (
            _json_sha256(list(training_sample_ids)) if training_sample_ids is not None else None
        ),
        "training_input_sha256": training_input_sha256,
        "configuration": {"seed": seed, "n_ctrl": n_ctrl, "n_bins": n_bins, "n_train_rows": int(len(rows))},
        "coverage": float(coverage),
        "missing": list(missing),
        "control_genes": control_map,
        "orthogonal_genes": list(orthogonal),
        "orthogonal_reference": orthogonal_reference,
        "citation": citation,
        "method": "control-gene module score; not the GSEA procedure in Saul et al. 2022",
    }
    return artifact


def validate_fit_artifact(artifact: dict) -> None:
    """Reject unsupported, corrupted, or internally inconsistent fit artifacts."""
    if not isinstance(artifact, dict) or artifact.get("format") != FIT_ARTIFACT_FORMAT:
        raise ScoreError("not a senescore fit artifact")
    if artifact.get("version") != FIT_ARTIFACT_VERSION:
        raise ScoreError("unsupported fit artifact version")
    signature = artifact.get("signature")
    features = artifact.get("feature_genes")
    controls = artifact.get("control_genes")
    config = artifact.get("configuration")
    orthogonal_genes = artifact.get("orthogonal_genes")
    orthogonal_reference = artifact.get("orthogonal_reference")
    if not isinstance(signature, list) or not signature or artifact.get("signature_sha256") != _json_sha256(signature):
        raise ScoreError("fit artifact signature hash is invalid")
    if not isinstance(features, list) or not features or artifact.get("feature_schema_sha256") != _json_sha256(features):
        raise ScoreError("fit artifact feature schema hash is invalid")
    if any(not isinstance(gene, str) or not gene for gene in signature + features):
        raise ScoreError("fit artifact gene symbols must be nonempty strings")
    if len(features) != len(set(features)) or len(signature) != len(set(signature)):
        raise ScoreError("fit artifact gene symbols must be unique")
    if (not isinstance(controls, dict) or not isinstance(config, dict)
            or not isinstance(orthogonal_genes, list) or not isinstance(orthogonal_reference, dict)):
        raise ScoreError("fit artifact is missing preprocessing fields")
    if any(not isinstance(gene, str) or not gene for gene in orthogonal_genes):
        raise ScoreError("fit artifact orthogonal gene symbols must be nonempty strings")
    if len(orthogonal_genes) != len(set(orthogonal_genes)):
        raise ScoreError("fit artifact orthogonal gene symbols must be unique")
    present = [gene for gene in signature if gene in set(features)]
    if set(controls) != set(present):
        raise ScoreError("fit artifact controls do not match the present signature genes")
    for gene, selected in controls.items():
        if not isinstance(selected, list) or not selected or any(item not in features for item in selected):
            raise ScoreError(f"fit artifact has invalid controls for {gene}")
    for gene, reference in orthogonal_reference.items():
        if gene not in features or gene not in orthogonal_genes or not isinstance(reference, list) or len(reference) != 2:
            raise ScoreError(f"fit artifact has invalid orthogonal reference for {gene}")
        if (any(isinstance(value, bool) or not isinstance(value, (int, float)) for value in reference)
                or not all(np.isfinite(value) for value in reference) or reference[1] <= 0):
            raise ScoreError(f"fit artifact has invalid orthogonal statistics for {gene}")
    if any(gene in features and gene not in orthogonal_reference for gene in orthogonal_genes):
        raise ScoreError("fit artifact is missing marker reference statistics")
    training_rows = artifact.get("training_rows")
    if (not isinstance(training_rows, list) or len(training_rows) != config.get("n_train_rows")
            or any(isinstance(row, bool) or not isinstance(row, int) or row < 0 for row in training_rows)
            or len(training_rows) != len(set(training_rows))):
        raise ScoreError("fit artifact training row metadata is invalid")
    if any(isinstance(config.get(key), bool) or not isinstance(config.get(key), int) or config[key] <= 0
           for key in ("n_ctrl", "n_bins", "n_train_rows")):
        raise ScoreError("fit artifact configuration is invalid")
    coverage = artifact.get("coverage")
    if isinstance(coverage, bool) or not isinstance(coverage, (int, float)) or not 0 < coverage <= 1:
        raise ScoreError("fit artifact coverage is invalid")
    if not isinstance(artifact.get("citation"), str) or not isinstance(artifact.get("method"), str):
        raise ScoreError("fit artifact citation or method metadata is invalid")
    training_ids = artifact.get("training_sample_ids")
    if training_ids is not None:
        if not isinstance(training_ids, list) or len(training_ids) != config["n_train_rows"]:
            raise ScoreError("fit artifact training sample IDs do not match the fitted row count")
        if any(not isinstance(sample, str) or not sample for sample in training_ids):
            raise ScoreError("fit artifact training sample IDs must be nonempty strings")
        if len(training_ids) != len(set(training_ids)):
            raise ScoreError("fit artifact training sample IDs must be unique")
        if artifact.get("training_sample_ids_sha256") != _json_sha256(training_ids):
            raise ScoreError("fit artifact training sample ID hash is invalid")
    elif artifact.get("training_sample_ids_sha256") is not None:
        raise ScoreError("fit artifact has a training sample ID hash without IDs")
    if artifact.get("training_input_sha256") is not None:
        digest = artifact["training_input_sha256"]
        if not isinstance(digest, str) or len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest.lower()):
            raise ScoreError("fit artifact input hash is invalid")


def transform_module_score(matrix: np.ndarray, gene_names: list[str], artifact: dict):
    """Apply a frozen fit artifact to a compatible later expression batch."""
    validate_fit_artifact(artifact)
    matrix, names, signature = _validated_names(
        matrix, gene_names, artifact["signature"], artifact["configuration"]["n_ctrl"],
        artifact["configuration"]["n_bins"],
    )
    if names != artifact["feature_genes"]:
        raise ScoreError("expression feature schema/order does not match the fit artifact")
    index = {gene: i for i, gene in enumerate(names)}
    present = [gene for gene in signature if gene in index]
    missing = [gene for gene in signature if gene not in index]
    sig_cols = np.array([index[gene] for gene in present])
    control_cols = [np.array([index[name] for name in artifact["control_genes"][gene]]) for gene in present]
    orthogonal_reference = {
        gene: tuple(values) for gene, values in artifact["orthogonal_reference"].items()
    }
    scores, orthogonal_dict = _apply_controls(
        matrix, names, present, sig_cols, control_cols,
        orthogonal_genes=artifact["orthogonal_genes"], orthogonal_reference=orthogonal_reference,
    )
    return {
        "scores": scores,
        "coverage": artifact["coverage"],
        "n_signature_present": len(present),
        "n_signature_missing": len(missing),
        "missing": missing,
        "configuration": dict(artifact["configuration"]),
        "control_genes": artifact["control_genes"],
        "orthogonal_z": orthogonal_dict,
        "orthogonal_fit_on": "fit_artifact_training_rows",
        "citation": artifact["citation"],
        "method": artifact["method"],
        "controls_fit_on": "fit_artifact_training_rows",
    }


def module_score(
    matrix: np.ndarray,
    gene_names: list[str],
    signature: tuple[str, ...] | list[str],
    seed: int = 0,
    n_ctrl: int = 5,
    n_bins: int = 20,
    block_from_controls: set[str] | None = None,
    orthogonal: tuple[str, ...] = ORTHOGONAL,
    citation: str = CITATION,
):
    matrix, names, signature = _validated_names(matrix, gene_names, signature, n_ctrl, n_bins)
    present, missing, coverage, sig_cols, control_cols = _select_controls(
        matrix.mean(axis=0), names, signature, seed=seed, n_ctrl=n_ctrl, n_bins=n_bins,
        block_from_controls=block_from_controls,
    )
    scores, orthogonal_dict = _apply_controls(matrix, names, present, sig_cols, control_cols, orthogonal_genes=orthogonal)
    return {
        "scores": scores,
        "coverage": coverage,
        "n_signature_present": len(present),
        "n_signature_missing": len(missing),
        "missing": missing,
        "configuration": {"seed": seed, "n_ctrl": n_ctrl, "n_bins": n_bins},
        "control_genes": {gene: [names[i] for i in cols] for gene, cols in zip(present, control_cols)},
        "orthogonal_z": orthogonal_dict,
        "orthogonal_fit_on": "all_rows",
        "citation": citation,
        "method": "control-gene module score; not the GSEA procedure in Saul et al. 2022",
    }


def module_score_train_only(
    matrix: np.ndarray,
    gene_names: list[str],
    signature: tuple[str, ...] | list[str],
    train_rows,
    seed: int = 0,
    n_ctrl: int = 5,
    n_bins: int = 20,
    block_from_controls: set[str] | None = None,
    orthogonal: tuple[str, ...] = ORTHOGONAL,
    citation: str = CITATION,
    *,
    gene_set_key: str | None = None,
    training_sample_ids: list[str] | None = None,
    training_input_sha256: str | None = None,
):
    """Fit bins and control genes on training rows only, then score every row.

    Cohort-wide means leak test-sample expression into the control set. This
    function freezes that choice on `train_rows`. It does not make the score a
    senescence assay.
    """
    rows = _validate_train_rows(train_rows, len(matrix))
    artifact = fit_module_score(
        matrix, gene_names, signature, rows, seed=seed, n_ctrl=n_ctrl, n_bins=n_bins,
        block_from_controls=block_from_controls, orthogonal=orthogonal, citation=citation,
        gene_set_key=gene_set_key, training_sample_ids=training_sample_ids,
        training_input_sha256=training_input_sha256,
    )
    result = transform_module_score(matrix, gene_names, artifact)
    result.update({
        "configuration": dict(artifact["configuration"]),
        "train_rows": list(artifact["training_rows"]),
        "orthogonal_fit_on": "train_rows_only",
        "controls_fit_on": "train_rows_only",
        "method": (
            "control-gene module score with bins and controls fit on training rows only; "
            "not the GSEA procedure in Saul et al. 2022"
        ),
        "fit_artifact": artifact,
    })
    return result

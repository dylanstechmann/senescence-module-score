"""Seurat-style control-gene module score, in numpy.

For each signature gene, subtract the average of control genes from the same
expression bin. A random gene set of the same size is the baseline, not a
second senescence clock.
"""

from __future__ import annotations

import numpy as np

from senescore.genes import CITATION, ORTHOGONAL, SENMAYO_HUMAN


class ScoreError(ValueError):
    pass


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


def pearson_correlation(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 2 or len(x) != len(y):
        return 0.0
    vx = x - x.mean()
    vy = y - y.mean()
    denom = float(np.sqrt(np.sum(vx ** 2) * np.sum(vy ** 2)))
    if denom == 0.0:
        return 0.0
    return float(np.sum(vx * vy) / denom)


def spearman_correlation(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 2 or len(x) != len(y):
        return 0.0
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
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


def _apply_controls(matrix, names, present, sig_cols, control_cols, *, z_reference_rows=None, orthogonal_genes=ORTHOGONAL):
    sig_mean = matrix[:, sig_cols].mean(axis=1)
    ctrl_mean = np.mean([matrix[:, cols].mean(axis=1) for cols in control_cols], axis=0)
    scores = sig_mean - ctrl_mean
    index = {gene: i for i, gene in enumerate(names)}
    orthogonal = {}
    for gene in orthogonal_genes:
        if gene in index:
            col = matrix[:, index[gene]]
            reference = col if z_reference_rows is None else col[z_reference_rows]
            std = float(reference.std()) or 1.0
            orthogonal[gene] = ((col - reference.mean()) / std).tolist()
    return scores, orthogonal


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
):
    """Fit bins and control genes on training rows only, then score every row.

    Cohort-wide means leak test-sample expression into the control set. This
    function freezes that choice on `train_rows`. It does not make the score a
    senescence assay.
    """
    matrix, names, signature = _validated_names(matrix, gene_names, signature, n_ctrl, n_bins)
    rows = np.asarray(train_rows)
    if rows.ndim != 1 or rows.dtype.kind not in "iu" or len(rows) < 2 or len(set(rows.tolist())) != len(rows):
        raise ScoreError("train_rows must be at least two unique sample indexes")
    if int(rows.min()) < 0 or int(rows.max()) >= len(matrix):
        raise ScoreError("train_rows out of range")
    present, missing, coverage, sig_cols, control_cols = _select_controls(
        matrix[rows].mean(axis=0), names, signature, seed=seed, n_ctrl=n_ctrl, n_bins=n_bins,
        block_from_controls=block_from_controls,
    )
    scores, orthogonal_dict = _apply_controls(
        matrix, names, present, sig_cols, control_cols, z_reference_rows=rows, orthogonal_genes=orthogonal
    )
    return {
        "scores": scores,
        "coverage": coverage,
        "n_signature_present": len(present),
        "n_signature_missing": len(missing),
        "missing": missing,
        "configuration": {"seed": seed, "n_ctrl": n_ctrl, "n_bins": n_bins, "n_train_rows": int(len(rows))},
        "train_rows": [int(i) for i in rows],
        "control_genes": {gene: [names[i] for i in cols] for gene, cols in zip(present, control_cols)},
        "orthogonal_z": orthogonal_dict,
        "orthogonal_fit_on": "train_rows_only",
        "citation": citation,
        "controls_fit_on": "train_rows_only",
        "method": (
            "control-gene module score with bins and controls fit on training rows only; "
            "not the GSEA procedure in Saul et al. 2022"
        ),
    }


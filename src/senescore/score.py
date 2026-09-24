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


def module_score(
    matrix: np.ndarray,
    gene_names: list[str],
    signature: tuple[str, ...] | list[str],
    seed: int = 0,
    n_ctrl: int = 5,
    n_bins: int = 20,
    block_from_controls: set[str] | None = None,
):
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
    index = {gene: i for i, gene in enumerate(names)}
    present = [gene for gene in signature if gene in index]
    missing = [gene for gene in signature if gene not in index]
    coverage = len(present) / len(signature)
    if coverage < 0.6:
        raise ScoreError(f"only {coverage:.0%} of the signature is in the table")
    means = matrix.mean(axis=0)
    order = np.lexsort((np.array(names), means))
    bins = np.empty(len(names), dtype=np.int64)
    # equal-count bins so a rare expression level still has controls
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
    sig_mean = matrix[:, sig_cols].mean(axis=1)
    ctrl_mean = np.mean([matrix[:, cols].mean(axis=1) for cols in control_cols], axis=0)
    scores = sig_mean - ctrl_mean
    orthogonal = {}
    for gene in ORTHOGONAL:
        if gene in index:
            col = matrix[:, index[gene]]
            std = float(col.std()) or 1.0
            orthogonal[gene] = ((col - col.mean()) / std).tolist()
    return {
        "scores": scores,
        "coverage": coverage,
        "n_signature_present": len(present),
        "n_signature_missing": len(missing),
        "missing": missing,
        "configuration": {"seed": seed, "n_ctrl": n_ctrl, "n_bins": n_bins},
        "control_genes": {gene: [names[i] for i in cols] for gene, cols in zip(present, control_cols)},
        "orthogonal_z": orthogonal,
        "citation": CITATION,
        "method": "control-gene module score; not the GSEA procedure in Saul et al. 2022",
    }


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


def random_signature(gene_names: list[str], size: int, seed: int, forbidden: set[str]) -> tuple[str, ...]:
    rng = np.random.default_rng(seed)
    pool = [gene for gene in gene_names if gene.upper() not in forbidden]
    pick = rng.choice(pool, size=size, replace=False)
    return tuple(str(gene) for gene in pick)

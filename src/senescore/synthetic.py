"""Spike SenMayo into a background transcriptome. For the bake-off only."""

from __future__ import annotations

import numpy as np

from senescore.genes import ORTHOGONAL, SENMAYO_HUMAN


def spiked_cohort(n_background: int = 1500, n_per_group: int = 30, shift: float = 1.2, seed: int = 0):
    rng = np.random.default_rng(seed)
    background = [f"BG{i}" for i in range(n_background)]
    genes = list(SENMAYO_HUMAN) + list(ORTHOGONAL) + background
    means = rng.uniform(0.5, 5.0, size=len(genes))
    n = n_per_group * 2
    matrix = means[None, :] + rng.normal(0.0, 0.35, size=(n, len(genes)))
    senescent = np.zeros(n, dtype=bool)
    senescent[n_per_group:] = True
    sig_idx = np.array([genes.index(gene) for gene in SENMAYO_HUMAN])
    matrix[np.ix_(senescent, sig_idx)] += shift
    labels = np.array(["reference"] * n_per_group + ["spiked"] * n_per_group)
    return matrix, genes, labels

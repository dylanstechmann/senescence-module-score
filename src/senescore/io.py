"""CSV: first row is gene symbols, first column is a sample id."""

from __future__ import annotations

import csv

import numpy as np

from senescore.score import ScoreError


def read_expression_csv(path: str):
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle))
    rows = [row for row in rows if row and any(cell.strip() for cell in row)]
    if len(rows) < 3:
        raise ScoreError("need a header and at least two samples")
    header = rows[0]
    genes = [gene.strip().upper() for gene in header[1:]]
    if not genes or any(not gene for gene in genes):
        raise ScoreError("missing gene symbols")
    if len(genes) != len(set(g.upper() for g in genes)):
        raise ScoreError("duplicate gene symbols")
    ids = []
    data = []
    for line, row in enumerate(rows[1:], 2):
        if len(row) != len(header):
            raise ScoreError(f"row {line}: ragged expression table")
        sample = row[0].strip()
        if not sample or sample in ids:
            raise ScoreError(f"row {line}: sample IDs must be nonempty and unique")
        ids.append(sample)
        try:
            data.append([float(value) for value in row[1:]])
        except ValueError as exc:
            raise ScoreError(f"row {line}: expression values must be numeric") from exc
    matrix = np.array(data, dtype=np.float64)
    if matrix.shape[1] != len(genes):
        raise ScoreError("ragged expression table")
    if not np.isfinite(matrix).all():
        raise ScoreError("expression values must be finite")
    return ids, genes, matrix

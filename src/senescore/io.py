"""CSV: first row is gene symbols, first column is a sample id."""

from __future__ import annotations

import csv

import numpy as np

from senescore.score import ScoreError


def read_expression_csv(path: str):
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle))
    if len(rows) < 3:
        raise ScoreError("need a header and at least two samples")
    header = rows[0]
    if header[0] == "":
        genes = header[1:]
    else:
        genes = header[1:]
    if len(genes) != len(set(g.upper() for g in genes)):
        raise ScoreError("duplicate gene symbols")
    ids = []
    data = []
    for row in rows[1:]:
        if not row or all(not cell.strip() for cell in row):
            continue
        ids.append(row[0])
        data.append([float(value) for value in row[1:]])
    matrix = np.array(data, dtype=np.float64)
    if matrix.shape[1] != len(genes):
        raise ScoreError("ragged expression table")
    return ids, genes, matrix

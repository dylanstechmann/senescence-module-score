"""CSV: first row is gene symbols, first column is a sample id."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from senescore.score import ScoreError


def read_expression_csv(path: str):
    data, digest = read_expression_snapshot(path)
    rows = list(csv.reader(data.decode("utf-8-sig").splitlines()))
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


def read_expression_snapshot(path: str) -> tuple[bytes, str]:
    """Read one immutable input snapshot and return its SHA-256 provenance."""
    import hashlib

    data = Path(path).read_bytes()
    return data, hashlib.sha256(data).hexdigest()


def read_expression_csv_with_hash(path: str):
    """Parse an expression CSV and bind the parse to the bytes that were read."""
    import hashlib

    data = Path(path).read_bytes()
    rows = list(csv.reader(data.decode("utf-8-sig").splitlines()))
    rows = [row for row in rows if row and any(cell.strip() for cell in row)]
    if len(rows) < 3:
        raise ScoreError("need a header and at least two samples")
    header = rows[0]
    genes = [gene.strip().upper() for gene in header[1:]]
    if not genes or any(not gene for gene in genes):
        raise ScoreError("missing gene symbols")
    if len(genes) != len(set(genes)):
        raise ScoreError("duplicate gene symbols")
    ids, data_rows = [], []
    for line, row in enumerate(rows[1:], 2):
        if len(row) != len(header):
            raise ScoreError(f"row {line}: ragged expression table")
        sample = row[0].strip()
        if not sample or sample in ids:
            raise ScoreError(f"row {line}: sample IDs must be nonempty and unique")
        ids.append(sample)
        try:
            data_rows.append([float(value) for value in row[1:]])
        except ValueError as exc:
            raise ScoreError(f"row {line}: expression values must be numeric") from exc
    matrix = np.asarray(data_rows, dtype=np.float64)
    if matrix.ndim != 2 or matrix.shape[1] != len(genes):
        raise ScoreError("ragged expression table")
    if not np.isfinite(matrix).all():
        raise ScoreError("expression values must be finite")
    return ids, genes, matrix, hashlib.sha256(data).hexdigest()

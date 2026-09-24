"""Score a CSV, or print the synthetic bake-off."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

from senescore.genes import SENMAYO_HUMAN
from senescore.io import read_expression_csv
from senescore.score import ScoreError, cohen_d, module_score, random_signature
from senescore.synthetic import spiked_cohort


def bakeoff(seed: int = 0) -> dict:
    matrix, genes, labels = spiked_cohort(seed=seed)
    real = module_score(matrix, genes, SENMAYO_HUMAN, seed=seed, block_from_controls=set(SENMAYO_HUMAN))
    decoy = random_signature(genes, len(SENMAYO_HUMAN), seed=seed + 1, forbidden=set(SENMAYO_HUMAN))
    base = module_score(matrix, genes, decoy, seed=seed, block_from_controls=set(SENMAYO_HUMAN))
    spiked = labels == "spiked"
    return {
        "senmayo_cohen_d": round(cohen_d(real["scores"][spiked], real["scores"][~spiked]), 3),
        "random_set_cohen_d": round(cohen_d(base["scores"][spiked], base["scores"][~spiked]), 3),
        "coverage": real["coverage"],
        "citation": real["citation"],
        "method": real["method"],
        "note": "The spike is synthetic. A large d here does not validate the set on your tissue.",
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="SenMayo control-gene module score")
    sub = parser.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    score = sub.add_parser("score")
    score.add_argument("csv")
    score.add_argument("--seed", type=int, default=0)
    score.add_argument("--controls", type=int, default=5)
    score.add_argument("--bins", type=int, default=20)
    args = parser.parse_args(argv)
    if args.cmd == "demo":
        json.dump(bakeoff(), sys.stdout, indent=2)
        sys.stdout.write("\n")
        return 0
    try:
        ids, genes, matrix = read_expression_csv(args.csv)
        result = module_score(matrix, genes, SENMAYO_HUMAN, seed=args.seed,
                              n_ctrl=args.controls, n_bins=args.bins)
    except (ScoreError, OSError) as exc:
        parser.error(str(exc))
    payload = {
        "input_sha256": hashlib.sha256(Path(args.csv).read_bytes()).hexdigest(),
        "configuration": result["configuration"],
        "control_genes": result["control_genes"],
        "samples": [
            {"id": sample_id, "senmayo_module_score": round(float(value), 4)}
            for sample_id, value in zip(ids, result["scores"])
        ],
        "coverage": result["coverage"],
        "missing": result["missing"],
        "orthogonal_z": {
            gene: [round(v, 3) for v in values]
            for gene, values in result["orthogonal_z"].items()
        },
        "citation": result["citation"],
        "method": result["method"],
        "not": "Not a biological-age clock, not a diagnosis, and not evidence that a compound is senolytic.",
    }
    json.dump(payload, sys.stdout, indent=2, allow_nan=False)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

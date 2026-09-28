from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

from senescore.genes import GENE_SETS, SENMAYO_HUMAN, get_gene_set
from senescore.io import read_expression_csv
from senescore.score import (
    ScoreError,
    cohen_d,
    module_score,
    module_score_train_only,
    pearson_correlation,
    random_signature,
    spearman_correlation,
)
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
    parser = argparse.ArgumentParser(description="Control-gene senescence module score")
    sub = parser.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")

    score = sub.add_parser("score", help="Score an expression CSV against a published gene set")
    score.add_argument("csv")
    score.add_argument("--gene-set", default="senmayo", choices=sorted(GENE_SETS.keys()),
                       help="published gene set to score (default: senmayo)")
    score.add_argument("--seed", type=int, default=0)
    score.add_argument("--controls", type=int, default=5)
    score.add_argument("--bins", type=int, default=20)
    score.add_argument("--train-samples", help="text file of sample ids; bins and controls are fit on these rows only")

    compare = sub.add_parser("compare", help="Compare two published gene sets head-to-head on the same expression table")
    compare.add_argument("csv")
    compare.add_argument("--set-a", default="senmayo", choices=sorted(GENE_SETS.keys()),
                         help="first gene set (default: senmayo)")
    compare.add_argument("--set-b", default="fridman", choices=sorted(GENE_SETS.keys()),
                         help="second gene set (default: fridman)")
    compare.add_argument("--seed", type=int, default=0)
    compare.add_argument("--controls", type=int, default=5)
    compare.add_argument("--bins", type=int, default=20)

    args = parser.parse_args(argv)

    if args.cmd == "demo":
        json.dump(bakeoff(), sys.stdout, indent=2)
        sys.stdout.write("\n")
        return 0

    if args.cmd == "compare":
        try:
            ids, genes, matrix = read_expression_csv(args.csv)
            set_a_info = get_gene_set(args.set_a)
            set_b_info = get_gene_set(args.set_b)
            res_a = module_score(
                matrix, genes, set_a_info["symbols"],
                seed=args.seed, n_ctrl=args.controls, n_bins=args.bins,
                block_from_controls=set(set_a_info["symbols"]),
                orthogonal=set_a_info["orthogonal"],
                citation=set_a_info["citation"],
            )
            res_b = module_score(
                matrix, genes, set_b_info["symbols"],
                seed=args.seed, n_ctrl=args.controls, n_bins=args.bins,
                block_from_controls=set(set_b_info["symbols"]),
                orthogonal=set_b_info["orthogonal"],
                citation=set_b_info["citation"],
            )
        except (ScoreError, OSError) as exc:
            parser.error(str(exc))

        p_corr = pearson_correlation(res_a["scores"], res_b["scores"])
        s_corr = spearman_correlation(res_a["scores"], res_b["scores"])

        payload = {
            "input_sha256": hashlib.sha256(Path(args.csv).read_bytes()).hexdigest(),
            "comparison": {
                "set_a": {
                    "key": args.set_a,
                    "name": set_a_info["name"],
                    "citation": set_a_info["citation"],
                    "n_signature_total": len(set_a_info["symbols"]),
                    "n_signature_present": res_a["n_signature_present"],
                    "coverage": round(res_a["coverage"], 4),
                    "mean_score": round(float(res_a["scores"].mean()), 4),
                    "std_score": round(float(res_a["scores"].std()), 4),
                },
                "set_b": {
                    "key": args.set_b,
                    "name": set_b_info["name"],
                    "citation": set_b_info["citation"],
                    "n_signature_total": len(set_b_info["symbols"]),
                    "n_signature_present": res_b["n_signature_present"],
                    "coverage": round(res_b["coverage"], 4),
                    "mean_score": round(float(res_b["scores"].mean()), 4),
                    "std_score": round(float(res_b["scores"].std()), 4),
                },
            },
            "correlation": {
                "pearson_r": round(p_corr, 4),
                "spearman_rho": round(s_corr, 4),
            },
            "samples": [
                {
                    "id": sample_id,
                    f"{args.set_a}_score": round(float(sa), 4),
                    f"{args.set_b}_score": round(float(sb), 4),
                    "diff": round(float(sa - sb), 4),
                }
                for sample_id, sa, sb in zip(ids, res_a["scores"], res_b["scores"])
            ],
            "configuration": {"seed": args.seed, "n_ctrl": args.controls, "n_bins": args.bins},
            "not": "Not a biological-age clock, not a diagnosis, and not evidence that a compound is senolytic.",
        }
        json.dump(payload, sys.stdout, indent=2, allow_nan=False)
        sys.stdout.write("\n")
        return 0

    try:
        ids, genes, matrix = read_expression_csv(args.csv)
        gs_info = get_gene_set(args.gene_set)
        if args.train_samples:
            wanted = [line.strip() for line in Path(args.train_samples).read_text(encoding="utf-8").splitlines() if line.strip()]
            missing = [sample for sample in wanted if sample not in ids]
            if missing or len(wanted) != len(set(wanted)):
                raise ScoreError("train sample ids must be unique and present in the expression table")
            lookup = {sample: i for i, sample in enumerate(ids)}
            result = module_score_train_only(
                matrix, genes, gs_info["symbols"], [lookup[sample] for sample in wanted],
                seed=args.seed, n_ctrl=args.controls, n_bins=args.bins,
                block_from_controls=set(gs_info["symbols"]),
                orthogonal=gs_info["orthogonal"],
                citation=gs_info["citation"],
            )
        else:
            result = module_score(
                matrix, genes, gs_info["symbols"], seed=args.seed,
                n_ctrl=args.controls, n_bins=args.bins,
                block_from_controls=set(gs_info["symbols"]),
                orthogonal=gs_info["orthogonal"],
                citation=gs_info["citation"],
            )
    except (ScoreError, OSError) as exc:
        parser.error(str(exc))
    payload = {
        "input_sha256": hashlib.sha256(Path(args.csv).read_bytes()).hexdigest(),
        "gene_set": args.gene_set,
        "configuration": result["configuration"],
        "control_genes": result["control_genes"],
        "samples": [
            {"id": sample_id, f"{args.gene_set}_module_score": round(float(value), 4)}
            for sample_id, value in zip(ids, result["scores"])
        ],
        "coverage": result["coverage"],
        "missing": result["missing"],
        "orthogonal_z": {
            gene: [round(v, 3) for v in values]
            for gene, values in result["orthogonal_z"].items()
        },
        "orthogonal_fit_on": result["orthogonal_fit_on"],
        "citation": result["citation"],
        "method": result["method"],
        "not": "Not a biological-age clock, not a diagnosis, and not evidence that a compound is senolytic.",
    }
    json.dump(payload, sys.stdout, indent=2, allow_nan=False)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


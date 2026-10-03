from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

from senescore.genes import GENE_SETS, SENMAYO_HUMAN, get_gene_set
from senescore.io import read_expression_csv_with_hash
from senescore.score import (
    ScoreError,
    cohen_d,
    module_score,
    module_score_train_only,
    pearson_correlation,
    random_signature,
    spearman_correlation,
    transform_module_score,
    validate_fit_artifact,
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

    score = sub.add_parser("score", help="Score an expression CSV against a gene panel")
    score.add_argument("csv")
    score.add_argument("--gene-set", default="senmayo", choices=sorted(GENE_SETS.keys()),
                       help="gene panel to score (default: senmayo)")
    score.add_argument("--seed", type=int, default=0)
    score.add_argument("--controls", type=int, default=5)
    score.add_argument("--bins", type=int, default=20)
    score.add_argument("--train-samples", help="text file of sample ids; bins and controls are fit on these rows only")
    score.add_argument("--fit-artifact", help="write frozen training controls and marker references as JSON")
    score.add_argument("--apply-artifact", help="apply a previously saved fit artifact to this expression CSV")

    compare = sub.add_parser("compare", help="Compare two gene panels on the same expression table")
    compare.add_argument("csv")
    compare.add_argument("--set-a", default="senmayo", choices=sorted(GENE_SETS.keys()),
                         help="first gene panel (default: senmayo)")
    compare.add_argument("--set-b", default="fridman", choices=sorted(GENE_SETS.keys()),
                         help="second gene panel (default: fridman)")
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
            ids, genes, matrix, input_sha256 = read_expression_csv_with_hash(args.csv)
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
            "input_sha256": input_sha256,
            "comparison": {
                "set_a": {
                    "key": args.set_a,
                    "name": set_a_info["name"],
                    "citation": set_a_info["citation"],
                    "source_status": set_a_info.get("source_status"),
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
                    "source_status": set_b_info.get("source_status"),
                    "n_signature_total": len(set_b_info["symbols"]),
                    "n_signature_present": res_b["n_signature_present"],
                    "coverage": round(res_b["coverage"], 4),
                    "mean_score": round(float(res_b["scores"].mean()), 4),
                    "std_score": round(float(res_b["scores"].std()), 4),
                },
            },
            "correlation": {
                "pearson_r": round(p_corr, 4) if p_corr is not None else None,
                "spearman_rho": round(s_corr, 4) if s_corr is not None else None,
                "undefined_reason": "non-finite, mismatched, or constant score vectors" if p_corr is None or s_corr is None else None,
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

    if args.apply_artifact and (args.train_samples or args.fit_artifact):
        parser.error("--apply-artifact cannot be combined with --train-samples or --fit-artifact")
    if args.fit_artifact and not args.train_samples:
        parser.error("--fit-artifact requires --train-samples")

    try:
        ids, genes, matrix, input_sha256 = read_expression_csv_with_hash(args.csv)
        if args.apply_artifact:
            artifact = json.loads(Path(args.apply_artifact).read_text(encoding="utf-8"))
            validate_fit_artifact(artifact)
            result = transform_module_score(matrix, genes, artifact)
            output_gene_set = artifact.get("gene_set_key") or "custom"
        else:
            gs_info = get_gene_set(args.gene_set)
            output_gene_set = args.gene_set
        if args.apply_artifact:
            pass
        elif args.train_samples:
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
                gene_set_key=args.gene_set,
                gene_set_source_status=gs_info.get("source_status"),
                training_sample_ids=wanted,
                training_input_sha256=input_sha256,
            )
            if args.fit_artifact:
                Path(args.fit_artifact).write_text(
                    json.dumps(result["fit_artifact"], indent=2, allow_nan=False) + "\n",
                    encoding="utf-8",
                )
        else:
            result = module_score(
                matrix, genes, gs_info["symbols"], seed=args.seed,
                n_ctrl=args.controls, n_bins=args.bins,
                block_from_controls=set(gs_info["symbols"]),
                orthogonal=gs_info["orthogonal"],
                citation=gs_info["citation"],
            )
    except (ScoreError, OSError, ValueError) as exc:
        parser.error(str(exc))
    payload = {
        "input_sha256": input_sha256,
        "gene_set": output_gene_set,
        "gene_set_source_status": (
            artifact.get("gene_set_source_status") if args.apply_artifact
            else gs_info.get("source_status")
        ),
        "configuration": result["configuration"],
        "control_genes": result["control_genes"],
        "samples": [
            {"id": sample_id, f"{output_gene_set}_module_score": round(float(value), 4)}
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
    if args.apply_artifact:
        payload["fit_artifact_version"] = artifact["version"]
        payload["training_sample_ids_sha256"] = artifact.get("training_sample_ids_sha256")
        payload["training_input_sha256"] = artifact.get("training_input_sha256")
    elif args.train_samples:
        payload["training_sample_ids_sha256"] = result["fit_artifact"]["training_sample_ids_sha256"]
        payload["fit_artifact_path"] = str(Path(args.fit_artifact)) if args.fit_artifact else None
    json.dump(payload, sys.stdout, indent=2, allow_nan=False)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

# Agent instructions — senescence-module-score

Work only in this repository. This scores the published SenMayo gene list
(Saul et al. 2022, 125 genes) with a control-bin average. It is not the paper's
GSEA. CDKN1A and CDKN2A are optional orthogonal z-scores because they are not
in SenMayo.

## Do not

- Call the score a biological age, a senolytic result, or a reason to take a compound.
- Drop the 60% gene-coverage refusal.
- Fit control bins on the full supervised table and describe that as training-only. Use `--train-samples` when the split matters.
- Change the gene list to “improve” the synthetic Cohen's d. The demo d near 30 is a software spike, not biology.

## First commands

```bash
python -m pip install -e .
python -m unittest discover -s tests -v
PYTHONPATH=src python3 -m senescore.cli demo
```

## Improve, in this order

1. If scoring math changes, extend tests for: duplicate symbols, NaN rejection, column reorder stability, and train-only controls not equal to full-cohort controls on a fixture where they should differ.
2. Optional: a one-page note on a **public** expression table you actually download outside git, with the input hash. Do not commit the matrix if it is large.
3. Do not rank molecules. Point compound questions at `geroscience-compound-atlas`.

## Done when

Demo still shows a large SenMayo d and a small random-set d, and tests pass.

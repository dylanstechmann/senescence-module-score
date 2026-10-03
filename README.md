# Senescence module score

A control-gene module score for the human **SenMayo** gene set (125 genes) published by Saul et al., Nature Communications 2022 ([10.1038/s41467-022-32552-1](https://doi.org/10.1038/s41467-022-32552-1)).

The symbols match the paper's supplementary table. The procedure does **not**. Saul et al. used GSEA. This repo uses a Seurat-style score: each signature gene is compared with control genes from the same expression bin, then averaged. CDKN1A and CDKN2A are reported as optional orthogonal z-scores when they are present, because they are **not** members of SenMayo.

## The bake-off

`senescore demo` spikes the 125 genes in half of a synthetic cohort and scores SenMayo against a size-matched random gene set. Controls for both sets are blocked out of the SenMayo genes, so the null is not punished just for accidentally using the spiked genes as controls.

On seed 0 that demo reports a SenMayo Cohen's d around 30 and a random-set d under 0.4. The unit test only requires d > 2 versus |random d| < 1. The huge synthetic d is what happens when 125 genes move together by a fixed offset. It is a software check. It is not an effect size in tissue, and it is not evidence that a compound is senolytic.

## Non-goals

- A biological-age clock
- A diagnosis, a senolytic recommendation, or a dose
- Training on private medical records
- Substituting for a senescence assay (SA-β-gal, p16 protein, whatever the lab actually validates)

For compounds and evidence grades, see [geroscience-compound-atlas](https://github.com/dylanstechmann/geroscience-compound-atlas). This repo scores a transcriptome table. It does not rank molecules.

## Run

```bash
make test
PYTHONPATH=src python3 -m senescore.cli demo
PYTHONPATH=src python3 -m senescore.cli score path/to/expression.csv --gene-set senmayo
PYTHONPATH=src python3 -m senescore.cli score path/to/expression.csv --gene-set fridman
PYTHONPATH=src python3 -m senescore.cli compare path/to/expression.csv --set-a senmayo --set-b fridman
```

CSV shape: a header row of gene symbols, and a sample id in column 1. Values should already be on a log-expression scale. If fewer than 60% of the chosen gene set is present, scoring refuses to invent the rest.

Supported gene sets:
- `senmayo`: Human SenMayo 125-gene signature (Saul et al. 2022)
- `fridman`: custom 66-gene panel inspired by senescence literature; not a source-transcribed Fridman signature
- `sasp`: custom 55-gene SASP-oriented panel; not a source-transcribed Coppé signature

Python 3.10+ and numpy.

## License

MIT. The gene set is the published SenMayo list; the paper remains the citation for the set.

## Traceable scoring (v0.2)

Install with `python -m pip install -e .`. Scoring now rejects duplicate/blank
sample IDs, duplicate normalized gene symbols, ragged tables, NaN and infinity.
Gene symbols are trimmed and uppercased consistently; the signature must be
nonempty and unique.

```bash
senescore score expression.csv --seed 0 --controls 5 --bins 20 > score.json
```

JSON includes the exact input hash, seed, bin/control settings and the actual
control genes sampled for each signature gene. With these controls the score
is reproducible and inspectable. Reordering input gene columns preserves
selection, including expression ties.

This is a cohort-dependent, per-signature-gene control average, not an exact
implementation of Seurat's pooled-control method. Bins and controls use means
across the supplied samples. Therefore scores from separately scored cohorts
are not automatically comparable, and scoring all samples before supervised
train/test splitting can leak cohort information. Raw counts are not normalized
by this command. Supply consistently processed log-expression values.

Cohen's d now reports an error for unequal constant groups, whose pooled
variance is zero, instead of incorrectly returning no effect. The synthetic
demo remains a software test, not a senescence diagnosis.

## Training-row controls (v0.3)

Bins and control genes depend on the samples you pass in. Scoring a full
supervised table before the split lets test rows choose the controls.

```bash
senescore score training.csv --train-samples train_ids.txt --fit-artifact senmayo-fit.json
senescore score followup_batch.csv --apply-artifact senmayo-fit.json
```

`train_ids.txt` is one sample id per line. Those rows alone set the expression
bins and the sampled control genes. Every row in the CSV is then scored with
that frozen set. Optional CDKN1A/CDKN2A z-scores also use the training rows
for centering and scaling; JSON records `orthogonal_fit_on`. Fitting
on all rows can leak held-out expression into those marker z-scores.
`--fit-artifact` writes the exact controls, marker reference statistics,
signature and feature-schema hashes, training sample IDs/hash, and source CSV
hash. `--apply-artifact` reuses those values without refitting; the later CSV
must have the same gene names in the same order or it fails closed.
`module_score_train_only`, `fit_module_score`, and `transform_module_score`
are the Python entry points. Fitting on every row still matches
`module_score`. A frozen artifact improves reproducibility; it does not make
scores comparable across independently collected cohorts or establish a
senescence diagnosis.
Training-row indexes must be distinct integers within the table; fractional
indexes are refused rather than silently truncated.

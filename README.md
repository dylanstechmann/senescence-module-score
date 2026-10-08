# Senescence module score

A control-gene module score for the human **SenMayo** gene set (125 genes) published by Saul et al., Nature Communications 2022 ([10.1038/s41467-022-32552-1](https://doi.org/10.1038/s41467-022-32552-1)).

This is a personal hobby and learning project, developed with substantial assistance from AI coding tools.

The symbols match the paper's supplementary table. The procedure does **not**. Saul et al. used GSEA. This repo uses a Seurat-style score: each signature gene is compared with control genes from the same expression bin, then averaged. CDKN1A and CDKN2A are reported as optional orthogonal z-scores when they are present, because they are **not** members of SenMayo.

## Real public-data check — GSE268487

The first external-data check found a **source-label ranking failure**: this project's
SenMayo module score ranked all three quiescent LF1 fibroblast libraries above
all three senescent libraries (directional AUROC 0.0; senescent-minus-quiescent
mean difference -0.413). Signed Fridman ranked these libraries in the expected
direction (AUROC 1.0), but 34 of 100 random panels also achieved AUROC 1.0.
Small within-condition variation makes perfect ranking easy in this dataset.

These are actual GEO expression counts, with bins and controls fitted only on
the three proliferating libraries. This post-2022 study provides an external
expression check; it has one cell line, no donor/experiment block metadata,
condition-specific processing differences and no linked functional assay per
RNA-seq library. It does not establish donor independence, tissue performance
or rejuvenation. The authors report inflammatory/SASP enrichment in their
nominally quiescent condition and question its suitability as a negative control;
source labels do not establish functional senescence-negative truth.
The reported ranking failure concerns this scoring procedure and
dataset, not a replication of the SenMayo paper's GSEA.

```bash
python -m pip install -e .
python -m senescore.public_validation --download
python -m unittest discover -s tests -v
```

The fixed plan, source receipts, individual scores, panel overlaps, marker
contrasts and interpretation are in [validation/GSE268487/REPORT.md](validation/GSE268487/REPORT.md).
The 634 KiB count file downloads into ignored `data/public/`; source and mapping
hashes are checked. Full fits, transformed expression and random panel details
are written into ignored `artifacts/`. The committed results include every
library's score and all contrasts; neither expression matrix is committed.

## Second prespecified check: GSE160356 (retinal endothelial cells)

The plan, pinned inputs and evaluator for a second external check were committed and pushed
([PLAN.md](validation/GSE160356/PLAN.md)) before any score existed. It uses nine HRMEC RNA-seq libraries (three
early passage, three late passage, three etoposide-treated), fits controls on the early-passage libraries only,
and scores each early-passage library with a fit from the other two so the baseline is not centred by its own fit.

Under the prespecified rule SenMayo **passed**: all six induced libraries scored above all three early-passage
libraries (directional AUROC 1.0 for each induction) and none of 100 random 125-gene panels matched its mean
difference. Read that narrowly. There are three libraries per condition (the smallest possible exact permutation p
is 0.1), passage is confounded with condition, and no clone or pairing is deposited. CDKN1A (p21) did not rise in
either induced condition. Post hoc, the late-passage shift is carried by about ten genes while several classic SASP
genes fell. Uniform random panels cannot show specificity against interferon or DNA-damage responses. The result
does not cancel the GSE268487 ranking failure above. The paper behind the data benchmarks its own signature against
earlier ones, so published results for established signatures may exist; this is not a blind test at the level of
the literature. Full report: [REPORT.md](validation/GSE160356/REPORT.md).

**Gene-removal follow-up (2026-10-08).** A frozen second plan removed 10 random SenMayo genes in each of 1,000 draws and re-applied the same verdict rule on the same nine libraries ([plan](validation/GSE160356/robustness/PLAN.md), [report](validation/GSE160356/robustness/REPORT.md)). The `pass` was kept in 956 draws and lost in 44 (4.4 %): in 43 the late-passage AUROC fell to 0.889, and in one the rank was 0.059. The etoposide contrast never lost it. So the earlier result is more fragile on the late-passage side than the headline suggests. It is the same libraries, not new data.

```bash
python -m senescore.hrmec_validation --download      # needs the pinned revision named in the report
python -m unittest discover -s tests -v
```

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
PYTHONPATH=src python3 -m senescore.cli compare path/to/expression.csv --set-a senmayo --set-b fridman_up --train-samples train_ids.txt
```

CSV shape: a header row of gene symbols, and a sample id in column 1. Values should already be on a log-expression scale. If fewer than 60% of the chosen gene set is present, scoring refuses to invent the rest.

Supported gene sets:
- `senmayo`: Human SenMayo 125-gene signature (Saul et al. 2022)
- `fridman`: legacy custom 66-gene panel; membership/direction are unverified
- `fridman_up`, `fridman_down`: source-pinned MSigDB C2 release 2025.1.Hs transcriptions of FRIDMAN_SENESCENCE_UP (77 genes, M9143) and FRIDMAN_SENESCENCE_DN (13 genes, M9487), both mapped to Fridman & Tainsky Table 2S. MSigDB content is CC BY 4.0; directional records include source-line hashes and are checked against the bundled two-line GMT snapshot. Their symbols and raw release snapshot match the Regen Workbench vendored copies.
- `fridman-signed`: compute the UP control-subtracted score minus the independently scored DOWN set. Use training-only fitting and apply the same frozen gene schema to later batches:

```bash
senescore fridman-signed expression.csv --train-samples train_ids.txt --fit-artifact fridman-fit.json
senescore fridman-signed followup.csv --apply-artifact fridman-fit.json
```

- `sasp`: custom 55-gene SASP-oriented panel; not a source-transcribed Coppé signature

The signed output is a directional expression score, not the original source study's analysis procedure and not an independently validated senescence assay. The Coppé study measured context-specific secreted proteins; this repository does not claim its custom SASP panel is a single canonical source signature.

Python 3.10+ and numpy.

## License

Software: MIT. SenMayo uses its published list and retains its paper citation.
The bundled MSigDB Fridman data use CC BY 4.0, as detailed in
[DATA_SOURCES.md](DATA_SOURCES.md).

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

## Audit follow-up — 2026-10-04

The Fridman UP/DN additions are retained and packaged with their pinned GMT
snapshot, including in wheels. The signed artifact now rejects mixing
directional components fitted on different training references, disagreement
with its wrapper, controls drawn from either directional set, inconsistent
coverage/missing-gene metadata, and duplicated controls. Each directional set
still independently requires at least 60% coverage.

Signed JSON carries the fitted reference IDs, schema hash, source metadata,
controls-fit basis and artifact hash. A fitting response hashes canonical JSON;
an apply response hashes the exact saved artifact bytes and labels that basis.
`compare --train-samples train_ids.txt` now fits both panels on those rows and
reports both control maps; its second-panel default is the source-pinned
`fridman_up` rather than the historical custom panel. Use an explicit
`--set-b fridman` to reproduce the old comparison. These are expression
contrasts and software checks; independent tissue/cell-type validation against
senescence and functional endpoints remains necessary.

Validation: **33 unittest cases** pass, including directional minimum coverage,
mixed-reference rejection and training-only comparisons. The installed-wheel
smoke check succeeds outside the source checkout. The seed-0 demo reports
SenMayo d=34.903 versus random-set d=0.294. These are synthetic software results.
See [data-source provenance and licensing](DATA_SOURCES.md).

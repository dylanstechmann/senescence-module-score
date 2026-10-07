# Prespecified external expression check: GSE160356 (HRMEC passage and etoposide)

**Verdict under the prespecified rule: `pass` — a narrow one.** SenMayo, scored with
controls fitted only on early-passage libraries, ranked all three late-passage libraries and all three
etoposide libraries above all three early-passage libraries (directional AUROC 1.0 for each contrast,
3 versus 3), and none of 100 uniform random 125-gene panels matched its mean difference (empirical rank
0.0099, the smallest possible with 100 panels). That is a descriptive ranking in one
retinal endothelial culture system with three libraries per condition, passage confounded with
condition, and no clone or pairing identifier. It is **not** an assay validation, not independence,
and not evidence of tissue performance. Three things temper it, all reported below: the late-passage
shift rests on about ten genes (post hoc), CDKN1A did not rise in either induced condition, and a
uniform random panel cannot show that the score is specific to senescence rather than to
inflammatory or DNA-damage responses. The earlier [GSE268487](../GSE268487/REPORT.md) check ended in a
ranking failure for the same score in lung fibroblasts; this result does not cancel it.

## What was frozen, and when

[PLAN.md](PLAN.md) was committed and pushed (senescence-module-score commit
`caa9eb38853449ffef7553eeda74f4509411ed8f`) **before** any score existed. A ResearchDesk plan freeze pinning that
commit, the split and the inputs was then recorded (`2026-10-07T19:54:19.269459+00:00`) and pushed in
regen-workbench commit `bb0a182`. The evaluator ran once from that clean committed revision
(`2026-10-07T19:54:38+00:00`, about five minutes of compute), with no re-run and no change to any panel, seed,
bin count, control count or threshold. No panel was refused for coverage. The
procedure matches the plan; see Deviations.

The freeze clock, the "results inspected" statement and the ledger actor are self-reported. The
independent evidence is server-side: GitHub Actions created its runs for the plan push at
2026-10-07T19:52:36Z and for the freeze-snapshot push at 19:54:34Z, before the evaluator's recorded start (19:54:38Z). The linked paper benchmarks its own signature
against earlier ones in these datasets, which I did not read before freezing, so published
results for established signatures may already exist.

## Source receipts

GEO [GSE160356](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160356): nine human
retinal microvascular endothelial cell (HRMEC) RNA-seq libraries, public 2023-10-27. Three early
passage (P3 to P8; GSM4872128, GSM4872131, GSM4872134), three late passage (P15 to P20; GSM4872129,
GSM4872132, GSM4872135) and three 1 µM etoposide for 4 days then 4 days of growth (P6 to P8;
GSM4872130, GSM4872133, GSM4872136). `GSE160356_RAW.tar`: 1,873,920 bytes, SHA-256
`93c1634dc2d6db494cf1592ceb14471e118acc1d99dd27dc1181f1bcea944770`; every member and the HGNC mapping are hashed in
[sources.json](sources.json). The series text says "three independent clones" compared with
replicative and etoposide-induced "counterparts", but no sample-level clone or pairing is deposited,
so none is used and all contrasts are unpaired. The linked paper
([Guduric-Fuchs et al. 2024](https://doi.org/10.1111/acel.14240), CC BY 4.0) reports cell-culture
senescence assays for these models, but no per-library functional value is deposited. GEO states no
dataset-specific license; the data are not relabelled MIT. `GSE160166` (ECFC array) was not used:
its platform table has no gene symbols.

The importer recognized 58,884 Ensembl rows (41,900 mapped to
symbols, 16,984 unmapped but kept in the CPM denominator); 5 HTSeq summary
rows were excluded and reported. SenMayo coverage was 0.992 (124 of 125 symbols).

## Primary result: SenMayo

Higher score = nominal senescence direction. Primary protocol: each early-passage library is scored
with a fit from the other two; late-passage and etoposide libraries use the fit from all three.
"Rank" is the two-sided empirical rank of the absolute mean difference among SenMayo and 100 uniform
random 125-gene panels (seeds 1000 to 1099). Exact permutation p over 20 assignments cannot be below 0.1.

| Contrast | Mean difference | Cohen's d | Directional AUROC | Exact permutation p | Rank vs random | Sensitivity: all-baseline fit (diff / AUROC) |
|---|---|---|---|---|---|---|
| Late passage − early passage | +0.177 | 2.82 | 1.000 | 0.1 | 0.0099 | +0.209 / 1.000 |
| Etoposide − early passage | +0.275 | 2.96 | 1.000 | 0.1 | 0.0099 | +0.308 / 1.000 |
| Etoposide − late passage (secondary) | +0.098 | 0.88 | 0.667 | 0.5 | 0.0099 | +0.098 / 0.667 |

| Panel / protocol | EP_CT_1 | EP_CT_2 | EP_CT_3 | LP_CT_1 | LP_CT_2 | LP_CT_3 | ETO_1 | ETO_2 | ETO_3 |
|---|---|---|---|---|---|---|---|---|---|
| SenMayo, primary | +0.225 | +0.207 | +0.206 | +0.386 | +0.303 | +0.479 | +0.474 | +0.364 | +0.625 |
| SenMayo, all-baseline fit | +0.162 | +0.192 | +0.186 | +0.386 | +0.303 | +0.479 | +0.474 | +0.364 | +0.625 |

Etoposide and late-passage libraries are not separated from each other (AUROC 0.667). The baseline
libraries score about +0.2, i.e. SenMayo genes average slightly higher than their expression-matched
controls even in early passage; only the contrasts are interpreted.

## Random panels

Under the same primary protocol, 100 uniform random 125-gene panels (excluding every registered panel
gene and CDKN1A/CDKN2A) gave, for late passage − early passage, mean differences from -0.104 to +0.086
(median -0.002) and a directional AUROC of 1.0 in
7 of 100; for etoposide − early passage, from -0.100 to +0.066
and AUROC 1.0 in 11 of 100. Perfect 3-versus-3 ranking happens by chance about 1 time in 20, so
AUROC alone is weak evidence here; SenMayo's mean differences (+0.177 and
+0.275) lie beyond every random panel. Uniform sets are not expression-matched or
function-matched, so this separates SenMayo from arbitrary gene sets, not from other biology such as
interferon, inflammatory or DNA-damage responses (the linked study itself reports an interferon signature in
these cells, and etoposide causes DNA damage).

## Other registered panels (primary protocol; mean difference / directional AUROC; descriptive only, no verdict)

| Panel | Late passage − early passage | Etoposide − early passage |
|---|---|---|
| Fridman UP (MSigDB, 77) | +0.135 / 0.889 | +0.275 / 1.000 |
| Fridman DOWN (MSigDB, 13; expected direction is lower) | -0.375 / 0.000 | -0.226 / 0.222 |
| Signed Fridman (UP − DOWN) | +0.475 / 1.000 | +0.533 / 1.000 |
| Fridman custom 66 (unverified) | +0.056 / 0.667 | +0.224 / 1.000 |
| SASP custom 55 (unverified) | +0.158 / 0.889 | +0.377 / 1.000 |

These panels overlap heavily, so agreement is not independent evidence: SenMayo shares 54 of the 55 SASP-custom
genes, 34 with the Fridman custom panel and 11 with Fridman UP (Pearson r across the nine libraries:
SenMayo with Fridman UP 0.96, with SASP custom 0.90). The DOWN set is expected to fall; it did in late passage
(AUROC 0.0) and, under the primary protocol, only partly for etoposide (AUROC 0.222; 0.0 under the all-baseline
fit); its leave-one-out baseline scores are widely spread (−0.218, −0.054, +0.262).

## Orthogonal expression markers

CDKN1A and CDKN2A are not SenMayo members. Reference mean and SD come from the three early-passage libraries.

| Marker | EP_CT_1 | EP_CT_2 | EP_CT_3 | LP_CT_1 | LP_CT_2 | LP_CT_3 | ETO_1 | ETO_2 | ETO_3 |
|---|---|---|---|---|---|---|---|---|---|
| CDKN1A log2(CPM+1) | 9.48 | 9.13 | 8.67 | 8.62 | 8.56 | 7.97 | 9.45 | 9.39 | 9.39 |
| CDKN1A z vs early passage | +1.2 | +0.1 | -1.3 | -1.4 | -1.6 | -3.4 | +1.1 | +0.9 | +0.9 |
| CDKN2A log2(CPM+1) | 4.37 | 3.36 | 3.87 | 4.23 | 3.54 | 5.49 | 4.94 | 3.78 | 4.67 |
| CDKN2A z vs early passage | +1.2 | -1.2 | +0.0 | +0.9 | -0.8 | +4.0 | +2.6 | -0.2 | +2.0 |

**CDKN1A (p21) did not rise in either induced condition**: it was lower in all three late-passage libraries (z −1.4,
−1.6, −3.4) and about level in the etoposide libraries (z +1.1, +0.9, +0.9). CDKN2A (p16) rose in some
induced libraries (LP_CT_3 z +4.0; ETO_1 +2.6, ETO_3 +2.0) but not in others (LP_CT_2 −0.8, ETO_2 −0.2). Neither marker
is a per-library functional assay and the series has one time point per library, so this
is a flag, not a refutation, but a uniform canonical-marker confirmation of senescence is not visible in these
libraries.

## Post hoc: which SenMayo genes move?

Written after the verdict and not part of it ([posthoc_gene_contributions.py](posthoc_gene_contributions.py),
[output](posthoc_gene_contributions.json)); it uses plain differences of log2(CPM+1) group means with no controls and
changes no panel. For late − early passage, 70 of 124 genes were higher, 39 lower and
the rest unchanged; the mean gene difference was +0.179 and **+0.005 without the ten largest
increases**. Those ten were PLAT (+3.35), IGFBP5 (+2.93), SERPINE2 (+2.76), MMP1 (+2.17), IGFBP6 (+2.13), IGFBP3 (+1.85), IL18 (+1.69), FAS (+1.63), ANGPTL4 (+1.59), CSF2RB (+1.51). Several classic SASP genes fell:
IGFBP2 (-4.21), CCL2 (-3.09), TNFRSF1B (-1.92), PGF (-1.01), ICAM1 (-1.01). For etoposide − early passage the shift is broader:
91 higher, 24 lower, mean +0.310 and +0.177 without the top ten
(MMP1 (+3.55), SERPINE2 (+2.18), IGFBP3 (+2.05), MMP10 (+2.01), FAS (+1.80), ANGPT1 (+1.49), ...). So the late-passage mean shift is carried by a small set of genes
(protease, protease-inhibitor and IGF-binding-protein genes feature among them) rather than by a general SASP program.
Whether the ranking itself would survive their removal was not tested.

## Procedural record

The receipt `evaluation_receipt.json` was bound to the earlier ResearchDesk freeze:
binding `bound_prospective`, every check passed (inputs pinned, split matches, no held-out library in
any fit, no overlap, method revision equal to the pinned commit, metric present, receipt created after the
freeze). The record's claim state is `exploratory_only` and it does **not** support a confirmatory
claim, because the freeze itself is `exploratory`: the grouping unit (a library) is not an identified
independent-unit level. The exported dossier verifies with scientific review `not_established_by_this_tool` and
reproduction `not_attempted`. It is in
[regen-workbench/studies/senescence-endothelial-challenge](https://github.com/dylanstechmann/regen-workbench/tree/main/studies/senescence-endothelial-challenge).

## Limits

- Three libraries per condition; the smallest possible exact permutation p is 0.1; no multiplicity correction.
- Passage differs by condition (and etoposide is a treatment, not a passage), no clone or pairing is deposited,
  processing batch is unknown, and there is no per-library functional assay.
- One retinal endothelial culture system. Nothing here concerns tissue, donors, patients, age or rejuvenation.
- The leave-one-out baseline protocol uses two libraries per fit, so baseline scores are noisy; the all-baseline
  fit is shown beside it and agrees in direction.
- The control-bin score is this repository's, not the GSEA procedure of the SenMayo paper.
- A pass on this rule means the rule was met, not that the score is validated.

## Deviations from the plan

None in the analysis. After the verdict I added the post hoc gene-level script and this report; neither
changes any prespecified quantity. The fixture-based tests, evaluator and plan were not edited after the freeze.

## Reproduce

```bash
git checkout caa9eb388534   # the frozen revision
python -m pip install -e .
python -m senescore.hrmec_validation --download
python validation/GSE160356/posthoc_gene_contributions.py
```

Generated 2026-10-07T19:54:38+00:00 by evaluator version 1 (numpy 1.26.4).
Raw archive and combined table stay in ignored `data/public/GSE160356/`; full fits and the transformed table in
ignored `artifacts/GSE160356/`.

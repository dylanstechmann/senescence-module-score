# Actual public-data expression check: GSE268487

**SenMayo failed the prespecified senescent-versus-quiescent ranking check.**
All three quiescent libraries received higher SenMayo module scores than all
three senescent libraries. Both conditions were higher than the proliferating
baseline. This fails the prespecified ranking of the source-labeled conditions;
it does not prove biological false positives. The authors report inflammatory/
SASP enrichment in the nominally quiescent condition and question its suitability
as a senescence control. Those labels are not independently verified
senescence-negative functional ground truth. This is a result for the
repository's control-subtracted score, not the SenMayo paper's GSEA procedure.

## Qualification and source receipts

Actual [GEO GSE268487](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268487)
contains 9 human LF1 lung fibroblast RNA-seq libraries: GSM8291859–61 are
proliferating, GSM8291862–64 quiescent, GSM8291865–67 replicatively senescent.
The study was deposited 2024-05-28 and publicly released 2025-07-31, after
SenMayo (2022) and Fridman/Tainsky (2008). The linked
[Dalgarno et al. preprint](https://pmc.ncbi.nlm.nih.gov/articles/PMC12262444/)
describes LF1 experiments. This qualifies as an external study expression
check of fixed published signatures. It does **not** prove donor/laboratory
independence from the literature that informed those signatures.

GSM titles supply reported replicate numbers, not donor or experiment blocks.
The publication explicitly calls its Hi-C replicates biological; this does
not establish that the RNA-seq libraries are independent donors. No donor IDs
or paired culture/experiment block IDs were deposited; they remain null.
Proliferating/senescent RNA was extracted with TRIzol, while quiescent pellets
were sent for processing, and the sequencing lane handling differed. These
processing differences confound condition. The linked publication is a
preprint. Its separate IMR90 RT-qPCR experiment is not the LF1 RNA-seq cohort
and was not joined as a functional endpoint. No linked numerical functional
senescence assay per RNA-seq library was available.

Raw count gzip: **649,147 bytes**, SHA-256
`ccf26d8b418f9c701c609195e6da0eb742accd5528e24d8811fb9f1fb6dc159a`.
The full SOFT metadata and HGNC downloads are also hashed in [sources.json](sources.json).
GEO provides public access without a dataset-specific license in this record;
the submitted data are not relabelled MIT. The linked article is CC BY-NC-ND
4.0. The compressed HGNC mapping facts are reusable under its
[no-restrictions data policy](https://www.genenames.org/about/), with attribution
to HUGO Gene Nomenclature Committee at the University of Cambridge.

## Frozen procedure and result

[PLAN.md](PLAN.md) was frozen before scores were calculated. The importer
recognized 58,780 Ensembl rows, mapping 41,895 locus rows to 41,859 symbols.
CPM denominators include unmapped Ensembl genes. HTSeq summaries, numeric/custom
identifiers and spike-ins are excluded and their totals reported. Counts sum
by symbol before `log2(CPM+1)`; measured zero genes remain. No differential
expression or signature optimization was used.

Bins, controls and marker references used **only proliferating GSM8291859–61**.
All 6 senescent/quiescent libraries were held out from control fitting.
Seed 0, 20 bins and 5 controls were fixed. No donor or experiment holdout is
claimed. The primary contrast below is **senescent minus quiescent**, 3 vs 3
held-out libraries; higher scores are ranked in the nominal senescence direction.

| Panel | Mapped coverage | Mean difference | Cohen's d | Directional AUROC | Exact two-sided permutation p |
|---|---:|---:|---:|---:|---:|
| SenMayo | 124/125 | -0.4125 | -15.731 | 0.000 | 0.10 |
| Fridman UP | 77/77 | +0.3657 | +7.893 | 1.000 | 0.10 |
| Fridman DOWN, uninverted | 13/13 | -0.5781 | -25.227 | 0.000 | 0.10 |
| Signed Fridman UP minus DOWN | 77/77; 13/13 | +0.8711 | +15.249 | 1.000 | 0.10 |
| Custom Fridman, unverified direction | 64/66 | -0.0959 | -0.980 | 0.111 | 0.20 |
| Custom SASP | 55/55 | -0.7971 | -24.172 | 0.000 | 0.10 |

No set was refused. SenMayo's missing symbol is CCL3L1; 111 of 124 measured
members had a nonzero count in at least one proliferating library, versus 116
in each held-out condition. Fridman UP had 76/77/76 expressed members in
proliferating/senescent/quiescent conditions; DOWN had 13/12/12. The custom
Fridman panel lacks CTGF and CYR61 under current approved-symbol mapping;
members and aliases were not changed to improve coverage.

Signed Fridman fits its UP and DOWN controls with their union blocked, so its
components differ from separately fitted directional panels; their exact
scores are included in [results.json](results.json). Mean SenMayo scores were
0.0311 proliferating, 1.0558 quiescent and 0.6433 senescent. Its separation from
proliferating libraries is a **calibration comparison**, not the held-out
source-label contrast.

With only 3 vs 3 libraries there are 20 label assignments; the smallest exact
two-sided p is 0.10. These unadjusted descriptive p-values assume exchangeability,
which processing confounding undermines. The large effect sizes reflect very
small within-condition variation and are not population effect estimates.

## Random panels, overlaps and expression markers

For 100 fixed 125-gene random panels (seeds 1000–1099), excluding all registered
panel members and CDKN1A/CDKN2A, mean differences ranged -0.1200 to +0.1143
(median -0.01885). SenMayo's absolute contrast exceeded all 100, giving an
empirical rank 1/101 = 0.00990, **but its direction was wrong**. Uniform random
panels are not expression matched and this rank does not establish specificity.
34 random panels achieved AUROC 1.0 and 52 achieved AUROC 0.0; perfect ranks
alone are weak evidence with these libraries. Every random panel's contrast is
saved in [random_scores.csv](random_scores.csv); reproducible members/individual
scores are retained in ignored artifacts.

SenMayo overlaps Fridman UP by 11 genes and the custom SASP set by 54 genes.
SenMayo/custom-SASP scores have Pearson r=0.990 across all 9 libraries; this
cannot be called independent assay agreement. All overlap memberships and
score correlations are in results.json.

CDKN1A and CDKN2A are excluded from SenMayo. Mean raw logCPM values were:

| Expression marker | Proliferating | Quiescent | Senescent | Held-out directional AUROC |
|---|---:|---:|---:|---:|
| CDKN1A | 8.1983 | 9.4236 | 10.0107 | 1.0 |
| CDKN2A | 1.5722 | 3.3263 | 5.5850 | 1.0 |

Their per-library logCPM and proliferating-reference z-scores are in results.json.
These are expression markers, not linked functional validation. Both CDKN1A and
CDKN2A overlap Fridman UP and the custom Fridman panel; neither independently
validates those panels. Both are orthogonal to the SenMayo membership here.
Independent review clarified this wording; the frozen plan and computed
overlapping_panel_keys were unchanged.

## Verification and reproduction

All **47 unittest cases passed** in Linux Docker, including the real-data import,
column/GSM matching, invalid integers, metadata conflicts, unmapped-gene library
denominators, duplicate loci, source hashes and baseline-only control fitting.
The original synthetic demo is unchanged (SenMayo d=34.903; random d=0.294).
That demo remains a software test.

Real-input verification exposed 45 valid pseudoautosomal-Y identifiers in the
count file. [IMPORTER_CORRECTION.md](IMPORTER_CORRECTION.md) records the correction
and first-run values without altering the frozen plan. All 45 PAR_Y rows have
zero counts here, so corrected library sizes and numerical scores are unchanged.
The final importer source hash and correction-receipt hash are in results.json.

```bash
python -m pip install -e .
python -m senescore.public_validation --download
python -m unittest discover -s tests -v
```

Run from the repository root. Source matrix downloads into ignored
`data/public/GSE268487/`; metadata/mapping hashes and the frozen plan hash are
verified. Full fit artifacts, transformed expression and random panel members
are written into ignored `artifacts/GSE268487/`. Only metadata facts, mapping
snapshot, source receipts and compact outcomes are tracked. A hash mismatch
fails rather than silently accepting an updated source or plan.

The useful next validation is new donors and another laboratory/experiment,
balanced library processing, quiescent and inflammatory controls, and functional
assays linked to the exact expression samples. These results do not demonstrate
biological age, senolysis, rejuvenation or a path to eternal youth.

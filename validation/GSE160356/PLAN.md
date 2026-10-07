# Frozen real-data evaluation plan: GSE160356 (HRMEC passage and etoposide)

Frozen on 2026-10-07, before any module score, contrast or random-panel result was
calculated on this series.

## What was and was not looked at before this freeze

Inspected: the GEO series and sample metadata (titles, characteristics, protocols,
supplementary file names), the series design text, file sizes and SHA-256 hashes, and the
identifier column of the nine count files (identical and in the same order in all nine).
Sixteen count lines from one early-passage library (the first eight and last eight rows of
GSM4872128, which include HTSeq summary rows) were displayed incidentally while checking the
file format. No other count value, score, contrast or differential-expression statistic was
calculated or viewed. A hash-and-structure check (`--validate-only`, which prints only row,
library and mapped-gene counts) may be run after this freeze and before scoring.

The linked article (Guduric-Fuchs et al., *Aging Cell* 2024, PMID 39422883,
doi:10.1111/acel.14240, CC BY 4.0) was read for design and qualification only: cell types,
treatments and which accessions it cites. Its figure captions state that it benchmarks its own
endothelial signature (EndoSEN) against earlier senescence signatures by GSEA in endothelial
datasets, which may include these libraries. I did not read any result for that benchmark.
**Published results for established signatures on this data may therefore exist, so this is not
a blind test at the level of the literature.** It tests this repository's control-bin score, as
specified below, and nothing else.

## Question

Do this repository's fixed, source-pinned gene panels, scored with controls fitted only on
early-passage libraries, rank late-passage and etoposide-treated human retinal microvascular
endothelial cell (HRMEC) libraries above early-passage libraries? If SenMayo does, is that
separation distinguishable from uniformly random 125-gene panels scored the same way? The
first GSE268487 check ended in a ranking failure for SenMayo in lung fibroblasts, so a failure
here is a live possibility and will be reported first if it happens.

## Qualification and limits

- GEO **GSE160356**, public 2023-10-27, submitted 2020-10-28: nine RNA-seq libraries,
  Illumina NovaSeq 6000, STAR and htseq-count against GRCh38 release 96. Three early passage
  (P3 to P8; GSM4872128, GSM4872131, GSM4872134), three late passage (P15 to P20; GSM4872129,
  GSM4872132, GSM4872135) and three 1 uM etoposide for 4 days then 4 days of growth (P6 to P8;
  GSM4872130, GSM4872133, GSM4872136). GSM titles give condition and replicate number.
- The series text says "three independent clones" compared with replicative and
  etoposide-induced "counterparts". No sample-level clone identifier or pairing is deposited.
  Replicate numbers must not be turned into clone or pair identifiers: clone, donor and
  experiment block stay null, and every contrast is unpaired.
- Passage range differs by condition, so passage is confounded with condition. The linked paper
  reports cell-culture-level senescence assays for these models, but no per-library functional
  value is deposited with the series, so none is joined.
- The accession postdates the SenMayo paper (2022) and Fridman and Tainsky (2008) as a public
  release; it does not establish independence from the literature that informed those panels.
  The linked study derived its own endothelial signature (EndoSEN) from related ECFC data and
  used these HRMEC libraries to apply it. EndoSEN is **not** scored here (it would be circular).
- Not used: GSE160166 (ECFC array). Its platform table (GPL21827) carries only 60-mer
  sequences and a genome build, with no gene symbols, so symbol-level scoring would need an
  alignment-based probe annotation that is not attempted here. This is an eligibility finding,
  not a result.

## Inputs and transformations

- The GEO `GSE160356_RAW.tar` (nine per-library `*.genename.htcounts.txt.gz` HTSeq tables),
  with the archive and every member pinned by SHA-256 in `sources.json`. Despite the file name,
  the first column is an Ensembl gene identifier.
- The pinned HGNC mapping snapshot already used for GSE268487
  (`../GSE268487/hgnc_mapping.tsv.gz`): approved, unambiguous, version-stripped Ensembl ID to
  current symbol; no aliases; no signature edits.
- Counts must be finite non-negative integers. HTSeq summary rows (`__no_feature`,
  `__ambiguous`, `__too_low_aQual`, `__not_aligned`, `__alignment_not_unique`) are excluded from
  the denominator and reported. The CPM denominator is the sum over **all** Ensembl gene rows
  including unmapped ones. Counts for one symbol are summed before `log2(CPM + 1)`. No expression
  or condition filter; measured zeros stay. This is a stated denominator choice, not a recovery
  of the study's own normalization.

## Fitting and scoring

Seed 0, 20 bins, 5 controls per gene, CDKN1A and CDKN2A blocked from ordinary controls, and
the 60 percent coverage refusal kept for every panel (a refusal is a result). Signatures are
unchanged: SenMayo (primary), source-pinned Fridman UP and DOWN, signed Fridman (UP minus
DOWN), and the explicitly custom Fridman and SASP panels (descriptive only).

Only the three early-passage libraries may select controls. Two protocols:

1. **Primary.** Each early-passage library is scored with a fit made from the *other two*
   (leave-one-library-out), so no baseline library helps choose the controls used to score it.
   Late-passage and etoposide libraries are scored with the fit from all three early-passage
   libraries. This matters because scoring the baseline with its own fit centres it near zero by
   construction.
2. **Sensitivity.** Every library is scored with the all-three early-passage fit. Reported
   beside the primary result, never substituted for it.

The held-out libraries never influence any fit; the evaluator counts and reports every library
a fit was allowed to see.

## Contrasts and the verdict, fixed in advance

Descriptive 3-versus-3 contrasts, higher score = nominal senescence direction: late passage
minus early passage and etoposide minus early passage (both **primary**), and etoposide minus late
passage (secondary). For each: mean difference, pooled Cohen's d, directional AUROC with ties,
and the exact two-sided label-permutation p over all 20 assignments (the smallest possible value
is 0.1, so none of this can reach conventional significance; no multiplicity correction; no
exchangeability claim given the passage confound).

Random check: 100 uniformly sampled 125-gene sets, seeds 1000 to 1099, from mapped genes
excluding the union of all registered panels and CDKN1A/CDKN2A (from signatures and controls),
scored under the primary protocol. The two-sided empirical rank of SenMayo's absolute mean
difference is (1 + number of random sets with an absolute difference at least as large) / 101,
for each primary contrast. Uniform sets are not expression-matched and do not establish
specificity.

Verdict, SenMayo only, primary protocol:

- `pass`: directional AUROC is 1.0 for both primary contrasts **and** the empirical rank is at
  most 0.05 for both.
- `ranking_only_not_distinguishable_from_random_panels`: both AUROCs are 1.0 but at least one
  rank is above 0.05.
- `fail`: any primary AUROC is below 1.0.
- `not_scored`: SenMayo was refused for low coverage.

No threshold, seed, panel member, bin count or control count may be changed after the first
scored run. Other panels are reported with the same statistics and carry no verdict.

## Reporting commitments

`REPORT.md` leads with the verdict, including a failure. It lists every library's score under
both protocols, all contrasts, refusals, panel overlaps, cross-score correlations, CDKN1A and
CDKN2A as orthogonal expression markers, the random-panel summary, and a Deviations section for
anything that differs from this plan. A pass would mean only that this one score ranks these
nine libraries in the stated direction in one retinal endothelial culture system with three
libraries per condition and a passage confound. It would not show tissue performance,
independence, causal senescence or rejuvenation, and it would not be an assay.

## Procedural record

The run emits `evaluation_receipt.json` in the generic
`regen-workbench/evaluation-receipt/1` form (producer revision, input hashes, the fold
`train = early passage, test = the six held-out libraries`, metric names, reported leakage) so
that a ResearchDesk plan freeze made before the run can be checked against it. That freeze
clock, the "results inspected" statement above and any ledger actor are self-reported.

## Reproducibility and publication

Raw archive and combined table stay in ignored `data/public/GSE160356/`; full fit artifacts and
the transformed table go to ignored `artifacts/GSE160356/`. Committed: this plan, the source
manifest with hashes, reduced sample metadata, the evaluator, compact results, per-seed random
scores, the receipt and the report. GEO states no dataset-specific license, so the submitted
data are not relabelled MIT. The linked article is CC BY 4.0; it is cited, not reproduced.
HGNC mapping facts are reused under the HGNC no-restrictions data policy with attribution to the
HUGO Gene Nomenclature Committee at the University of Cambridge.

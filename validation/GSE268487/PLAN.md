# Frozen real-data evaluation plan

Frozen on 2026-10-04 before calculating any module score. Source qualification,
count-table header/identifier inspection and metadata review preceded this plan.
Neither differential-expression statistics nor score results guided these choices.

## Qualification and limits

Use GEO **GSE268487**, public 2025-07-31, deposited 2024-05-28. It contains
three RNA-seq libraries each from proliferating, quiescent and replicatively
senescent human LF1 lung fibroblasts. GSM titles supply condition and replicate
numbers. The associated 2025 preprint describes new LF1 experiments; this
accession postdates Saul et al. 2022 (SenMayo) and Fridman/Tainsky 2008.
That provides an external study check of the fixed signatures. It does not
prove donor/laboratory independence: LF1 was used in older senescence research,
and the deposited metadata do not identify donors or independent experiment
blocks. Replicate numbers must not be invented into donor or paired block IDs.

Proliferating and senescent RNA was extracted using TRIzol; quiescent cell
pellets were sent for processing. Proliferating/senescent libraries spanned two
HiSeq lanes; quiescent libraries used one lane. These differences confound
condition. The study is a preprint. No linked numerical functional senescence
assay per RNA-seq library was identified. The separate IMR90 RT-qPCR experiment
must not be joined to LF1 RNA-seq samples.

## Inputs and transformations

- Actual GEO `GSE268487_pro_qui_sen_ft_counts.txt.gz`, not DE output or simulation.
- Actual GEO SOFT metadata; match count columns to GSM title/replicate exactly.
- HGNC complete set snapshot: use only approved, unambiguous Ensembl-to-symbol
  mappings. Strip Ensembl version suffixes; do not alias or modify signatures.
- Validate all counts as finite nonnegative integers. Exclude HTSeq summary
  rows, numeric/custom identifiers and spike-ins from the Ensembl gene library
  denominator; report their counts separately. This is an explicit denominator
  choice, not a claim to recover the study's edgeR normalization.
- CPM denominator includes **all** Ensembl gene rows, including unmapped genes.
  Mapped counts for the same symbol are summed before `log2(CPM + 1)`.
  Do not filter on expression or condition differences; retain measured zeros.
  Report mapped coverage and numbers expressed in each condition separately.

## Fitting and contrasts

Fit bins, sampled controls and CDKN1A/CDKN2A reference mean/SD using the three
**proliferating libraries only**. All quiescent and senescent libraries are
held out from control fitting. Freeze seed=0, bins=20, controls=5. Block
CDKN1A/CDKN2A from ordinary panel controls. Signed Fridman uses its existing
UP/DOWN union blocking behavior.
This is unsupervised baseline calibration, not a trained diagnostic classifier.
No donor/experiment holdout can be claimed for this single-line dataset.

Primary descriptive contrast: senescent minus quiescent (3 vs 3 held-out
libraries). Report individual scores, mean difference, pooled Cohen's d,
directional AUROC with ties, and exact two-sided label-permutation p over all
20 assignments. Small sample count and condition/processing confounding limit
inference. No threshold optimization, training accuracy, or diagnostic claims.
Senescent-minus-proliferating and quiescent-minus-proliferating contrasts are
secondary descriptive comparisons against the fitted baseline, not independent
test contrasts. No replicate-number pairing or donor bootstrap.

Score unchanged SenMayo, source-pinned Fridman UP and DOWN, signed UP-minus-DOWN,
and the existing explicitly custom Fridman and SASP panels. Keep the 60%
coverage refusal separately for each component; record a refusal as a result.
Report overlaps and cross-score correlations rather than independent-assay
agreement when panels share genes. CDKN1A/CDKN2A are orthogonal only to SenMayo;
CDKN1A overlaps Fridman. Report raw logCPM and baseline z-scores, without calling
these functional assay validation.

Random check: 100 uniformly sampled 125-gene sets from mapped genes, seeds
1000..1099; exclude the union of all registered panel genes and CDKN1A/CDKN2A
from signatures and controls. Compare signed mean differences and AUROC;
report the two-sided empirical rank of the absolute SenMayo mean difference.
Uniform sets are not expression-matched and do not establish specificity.
No seed selection or panel revision following this run.

## Reproducibility and publication

Raw matrix/full SOFT/HGNC snapshot remain in ignored `data/public/GSE268487/`.
Commit the prespecified plan, source URLs/hash manifest, reduced sample metadata,
an unambiguous HGNC mapping snapshot, importer/evaluator and compact results.
Store full fit artifacts and transformed table in ignored `artifacts/`.
Verification must test integer rejection, identifier ambiguity/duplicates,
denominator handling, correct GSM/column matching, actual training IDs,
held-out expression changes leaving fitted controls unchanged, and reproducibility.

GEO data: publicly accessible, no dataset-specific license asserted in the GEO
record. Do not relabel submitted data as MIT. HGNC permits reuse with requested
attribution to HUGO Gene Nomenclature Committee at the University of Cambridge.
The associated preprint is CC BY-NC-ND 4.0; link/cite it rather than redistributing
its text or figures. Source snapshots contain metadata facts, not article text.

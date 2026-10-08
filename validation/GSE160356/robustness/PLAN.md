# Frozen robustness plan: SenMayo gene removal on GSE160356

Frozen on 2026-10-08, before any reduced-panel score was calculated. This is a follow-up to the
GSE160356 check in `../PLAN.md` and `../REPORT.md` (verdict `pass`, narrowly: empirical rank 0.0099
for all three contrasts, the lowest value 100 random panels allow). That report also noted two
caveats this plan tests directly: CDKN1A did not rise, and a post-hoc look suggested a small number
of genes carry the late-passage shift.

This plan was written by an AI coding assistant. The owner has not reviewed it.

## What was and was not looked at before this freeze

The earlier run's results, report and post-hoc gene table were read; they are the reason for this
question. No reduced-panel score, random-removal result or any statistic from this plan has been
calculated. The data are the same nine libraries as before, so this is **not** a new test of the
score on independent data. It asks only how much the earlier result depends on which genes are in
the panel.

A timing check ran the zero-removal control (gap 0.0 against the committed scores) and the scoring path
for draw seed 2000 once before this freeze. For that draw only the coverage (1.0 of the supplied
symbols) and the elapsed time (about 2.4 s) were displayed; no score, contrast or rank was viewed.

## Question

If 10 randomly chosen SenMayo genes are removed, how often does the earlier verdict survive, and
how large is the late-passage and etoposide shift in the worst draws?

## Design, fixed in advance

- Input: `artifacts/GSE160356/expression_logCPM.csv` (the transformed table the earlier run wrote;
  ignored by git), checked first by a zero-removal control (below).
- Panel: SenMayo symbols present in the table (124 of 125 in the earlier run). Per draw, remove
  exactly 10 of them, uniformly without replacement, using `numpy.random.default_rng(seed)` with
  draw seeds **2000 to 2999** (1,000 draws). Remaining genes are the panel (114 genes).
- Scoring: identical to the earlier primary protocol: seed 0, 20 bins, 5 controls per gene,
  CDKN1A and CDKN2A blocked from controls, leave-one-library-out for early-passage libraries and
  the all-three fit for held-out libraries. The ten removed genes are **also kept out of the
  control pool**, so the only change is panel membership.
- Control: with nothing removed, the scores must equal the committed SenMayo primary scores in
  `../results.json` to within 1e-9. If not, the run stops and reports that.
- Per draw: mean difference and directional AUROC for the two primary contrasts, and the
  two-sided empirical rank of the absolute mean difference against the 100 committed
  random-panel values in `../random_scores.csv` (125-gene panels; the reduced panel has 114 genes,
  so the comparison is approximate and said to be).
- Per-draw verdict: the earlier rule, unchanged (`hrmec_validation.verdict`): `pass` needs both
  primary AUROCs equal to 1.0 and both ranks at most 0.05.

## What will be reported, whatever it shows

Number of draws that keep `pass`; number with both primary AUROCs equal to 1.0; the
distribution (min, 2.5 %, median, 97.5 %, max) of each primary mean difference; the draws with
the smallest mean differences and which genes they removed; and the highest rank observed. The
report leads with the draw count that loses `pass`, including if that is most of them.

No threshold, seed, removal size, bin count, control count or rule is changed after the first run.
The run happens once, from a clean committed revision, and `results.json` records that revision.

## Limits

- Same nine libraries, one retinal endothelial culture system, three libraries per condition and
  a passage confound. Surviving removal shows robustness to panel membership in these data only.
- Random removal does not test whether a *particular* set of genes drives the result; the earlier
  post-hoc table already points at candidates, and a targeted removal would be a different,
  separately frozen question.
- No Desk freeze record is made for this plan; the commit time and hash of this file are the
  record, and "results inspected" above is self-reported.
- No tissue generalization, rejuvenation, diagnostic or biological-age claim.

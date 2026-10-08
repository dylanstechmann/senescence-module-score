# Frozen targeted-removal plan: SenMayo on GSE160356

Frozen on 2026-10-08 before any targeted-removal score was calculated. Written by an AI coding
assistant; the owner has not reviewed it. Follows `../robustness/REPORT.md`, which found the
late-passage result fragile under random removal and listed PLAT, CSF2RB and IL6 as recurring
genes in losing draws (post hoc).

## What was looked at before this freeze

The random-removal results, the post-hoc recurrence table and the earlier post-hoc table of the
ten highest-moving genes per contrast (`../posthoc_gene_contributions.json`) were read. They are
where the gene lists below come from, so **these removals are chosen after seeing the data and
test dependence on those genes, not a new validation.** No targeted-removal score has been
calculated.

## Four fixed removals (scored once each, same protocol as the earlier primary run)

- `recurring3`: PLAT, CSF2RB, IL6.
- `late_top10`: the ten genes with the largest late-passage minus early-passage log2(CPM+1)
  difference in the earlier post-hoc table: PLAT, IGFBP5, SERPINE2, MMP1, IGFBP6, IGFBP3, IL18,
  FAS, ANGPTL4, CSF2RB.
- `etoposide_top10`: MMP1, SERPINE2, IGFBP3, MMP10, FAS, ANGPT1, PLAT, ANGPTL4, CXCL8, CSF2RB.
- `union_top`: the union of the two top-10 lists and recurring3.

Removed genes also stay out of the control pool. Scoring settings, the verdict rule and the 100
committed random-panel deltas are unchanged. The reduced panels are smaller than the 125-gene
random panels, so ranks are approximate.

## Reported whatever it shows

Per removal: panel size, both primary mean differences and AUROCs, ranks, and the verdict. The
report leads with how many of the four lose `pass`. No other removal is added after the run;
a different list needs a new frozen plan.

## Limits

Same nine libraries. Because the genes were picked from these same data, a loss shows how much
the result leans on them; a pass does not show the score is independent of them.

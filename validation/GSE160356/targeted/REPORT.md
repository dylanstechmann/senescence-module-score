# Targeted-removal report: SenMayo on GSE160356 (run 2026-10-08)

Personal hobby and learning project with substantial AI assistance; plan and run by an AI coding
assistant, not reviewed by the owner. Plan: [PLAN.md](PLAN.md), frozen in `e82d596`; the run
recorded that revision with no tracked-tree changes. The gene lists were chosen after seeing
these data (post-hoc table and the random-removal recurrence), so this measures dependence on
those genes. It is not a validation.

## Result

**Three of the four fixed removals lose the earlier `pass`.** The late-passage result depends
on a small set of genes; the etoposide result does not.

| Removal | Panel | Late minus early (AUROC, rank) | Etoposide minus early (AUROC, rank) | Verdict |
|---|---:|---|---|---|
| none (earlier run) | 124 | 0.177 (1.0, 0.0099) | 0.275 (1.0, 0.0099) | pass |
| `recurring3` PLAT, CSF2RB, IL6 | 121 | 0.124 (1.0, 0.0099) | 0.204 (1.0, 0.0099) | pass |
| `late_top10` | 114 | **0.014 (0.556, 0.72)** | 0.173 (1.0, 0.0099) | fail |
| `etoposide_top10` | 114 | **0.040 (0.778, 0.32)** | 0.122 (1.0, 0.0099) | fail |
| `union_top` (14 genes) | 110 | **-0.036 (0.333, 0.35)** | 0.108 (1.0, 0.0099) | fail |

Removing the ten genes with the largest late-passage shift leaves almost no late-passage
separation (mean difference 0.014, AUROC 0.556 against 1.0), and the union of the lists reverses
its sign. Those ten genes overlap heavily with the etoposide list (PLAT, SERPINE2, MMP1, IGFBP3,
FAS, ANGPTL4, CSF2RB), which is why both lists hurt the late-passage contrast. The three
recurring genes alone are not enough to lose the verdict.

## Reading

The earlier `pass` says SenMayo ranks these nine libraries in the nominal direction. This run
shows that, for late passage, the ranking is carried by about ten genes, several of them
extracellular-matrix and secreted factors (MMP1, SERPINE2, IGFBP3 and others) that many stress
and culture conditions would move. The score on this series should not be read as a broad
senescence readout. The etoposide contrast stays at AUROC 1.0 with rank 0.0099 after every
removal, so it is spread over more genes. Same nine libraries; no tissue, rejuvenation or
diagnostic claim.

## Deviations

None. Ranks compare smaller panels against the 125-gene random panels and are approximate.

# Robustness report: SenMayo gene removal on GSE160356 (run 2026-10-08)

This is a personal hobby and learning project, developed with substantial assistance from AI
coding tools. The plan was written and the run executed by an AI coding assistant; the owner has
not reviewed either. Plan: [PLAN.md](PLAN.md), frozen in commit `52a4604` before any draw was
scored. Code revision recorded in `results.json`: `52a4604d636360a98bece895e91296ae0c74f2c8`,
no tracked-tree changes.

## Result

**The earlier `pass` does not survive every panel. Removing 10 random SenMayo genes lost the
verdict in 44 of 1,000 draws (4.4 %); 956 draws kept `pass`.**

| Verdict across 1,000 draws | Draws |
|---|---:|
| `pass` | 956 |
| `fail` (a primary AUROC below 1.0) | 43 |
| `ranking_only_not_distinguishable_from_random_panels` | 1 |

- All 43 `fail` draws have the same cause: the **late-passage versus early-passage** AUROC fell
  from 1.0 to 0.889, which is one inverted library pair out of nine. The etoposide contrast kept
  an AUROC of 1.0 in all 1,000 draws and a rank of 0.0099 (the lowest 100 random panels allow) in
  every draw.
- The single `ranking_only` draw had both AUROCs at 1.0 and a late-passage rank of 0.059, just
  above the 0.05 rule.
- Mean differences (log2 CPM+1 scale, same transformation as before). Late passage minus early:
  full panel 0.177; across draws min 0.072, 2.5 % 0.111, median 0.174, 97.5 % 0.235, max 0.266.
  Etoposide minus early: full panel 0.275; min 0.161, 2.5 % 0.200, median 0.255, 97.5 % 0.303,
  max 0.336. Both contrasts keep their direction in every draw.
- Zero-removal control: the scoring path reproduced the committed SenMayo primary scores with a
  maximum absolute gap of 0.0.

## What this means, and does not

The etoposide result is robust to losing ten random genes. The late-passage result is more
fragile: about one panel in twenty-three (43 of 1,000, plus the near miss) loses a perfect
ranking, and the late-passage shift can fall to less than half its full-panel size. That fits the
earlier post-hoc observation that a modest number of genes carry the late-passage difference, and
it is a weaker statement than the original verdict made. It is the same nine libraries, so this
does not test the score on new data. Random removal also does not say which genes matter.

## Post-hoc look (not part of the frozen plan)

`posthoc_removed_genes.py` (output in `posthoc_removed_genes.json`) asks which removed genes
recur in the 44 losing draws. Against an overall loss rate of 4.4 %, draws that removed **PLAT**
lost `pass` in 18 of 109 (16.5 %), **CSF2RB** in 12 of 78 (15.4 %) and **IL6** in 10 of 76
(13.2 %). These are counts from overlapping random draws, not independent tests, and no
correction is made. They are leads for a separately frozen targeted-removal plan, nothing more.

## Deviations from the plan

None in the design, seeds, scoring settings or verdict rule. Disclosures: a one-draw timing check
(seed 2000) and the zero-removal control ran before the freeze, and only coverage and elapsed time
were displayed for that draw (stated in the plan). No Desk freeze record was made for this plan;
the commit is the record. The reduced panels have 114 genes while the 100 random comparison
panels have 125, so the ranks are approximate, as the plan said.

## Files

`PLAN.md`, `results.json` (summary, hashes, control), `draws.csv` (every draw, its removed genes,
both primary contrasts, ranks and verdict), `posthoc_removed_genes.py` and `.json`, this report.
Reproduce: `PYTHONPATH=src python3 -m senescore.hrmec_robustness` from the repository root, with
`artifacts/GSE160356/expression_logCPM.csv` from the earlier run (about 40 minutes on one core).
The frozen run writes once and refuses to run again where `results.json` exists.

## Limits

Same nine libraries and one retinal endothelial culture system, three libraries per condition and
a passage confound. No tissue generalization, rejuvenation, diagnostic or biological-age claim.

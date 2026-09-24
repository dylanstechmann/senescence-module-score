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
PYTHONPATH=src python3 -m senescore.cli score path/to/expression.csv
```

CSV shape: a header row of gene symbols, and a sample id in column 1. Values should already be on a log-expression scale. If fewer than 60% of SenMayo is present, scoring refuses to invent the rest.

Python 3.10+ and numpy.

## License

MIT. The gene set is the published SenMayo list; the paper remains the citation for the set.

# Importer correction discovered by real-input verification

The first scoring execution recognized only `ENSG` plus digits and an optional
version. It incorrectly excluded 45 valid `..._PAR_Y` annotation rows as custom
identifiers. A subsequent real-input integration test with strict malformed-ID
refusal exposed `ENSG00000002586.19_PAR_Y`. The frozen plan already required
**all Ensembl gene rows** in the library denominator and symbol-level summation;
the implementation was corrected to honor that plan.

X and pseudoautosomal-Y loci remain distinct count identifiers, while their
mapped counts sum into the same HGNC symbol before log transformation. Duplicate
versions of the same locus still fail. No signature, coverage rule, split,
control seed, bin count, random seeds or contrast was changed.

The first (incorrect-import) result had SenMayo senescent-minus-quiescent
delta -0.412502, directional AUROC 0.0, and signed Fridman delta +0.871142,
directional AUROC 1.0. Those values were inspected before the correction;
the correction is a source-format bug fix, not a result-selected analysis.
The final results and report use the corrected importer and record its version.

## Metadata retention correction

Independent review found that the extraction-protocol whitelist used
`!Sample_extraction_protocol_ch1` rather than GEO's actual
`!Sample_extract_protocol_ch1`. Version 3 retains the exact two values for every
sample, including the extraction, library kit and condition-specific lane
handling. The factual sample snapshot and its hash were regenerated from the
same full SOFT input. A regression test checks exact protocol retention, and
the real-input integration test checks the snapshot against the full source.
This changes metadata receipts only; the frozen plan, expression input,
signature membership, splits, seeds and numerical scores are preserved.

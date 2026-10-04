# Gene-set data provenance

The new Fridman sets are pinned MSigDB C2 release 2025.1.Hs transcriptions:
[`FRIDMAN_SENESCENCE_UP`](https://www.gsea-msigdb.org/gsea/msigdb/cards/FRIDMAN_SENESCENCE_UP)
(M9143, 77 members) and
[`FRIDMAN_SENESCENCE_DN`](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/FRIDMAN_SENESCENCE_DN.html)
(M9487, 13 members). The source cards were checked on October 4, 2026.
Both identify human NCBI Gene IDs and source Table 2S of
[Fridman and Tainsky, Oncogene 2008, PMID 18711403](https://pubmed.ncbi.nlm.nih.gov/18711403/).
This is a review of expression profiling studies, not a validated diagnostic.

The JSON fixtures retain each published database direction, systematic ID,
release, source table, mapped symbols and aligned Entrez IDs. Symbols are the
MSigDB release mappings; this project does not apply an additional alias map.
The bundled two-line GMT contains the release's directional symbol entries.
Import verifies exact ordered GMT/JSON membership and each raw-line SHA-256.
The normalized-membership hash is SHA-256 of canonical JSON containing the
ordered `symbols` and `entrez_ids` arrays (sorted keys, compact separators,
ASCII encoding). It binds the symbol-to-identifier mapping as well as membership.

These MSigDB gene sets are licensed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
Attribution: Broad Institute, Inc., Massachusetts Institute of Technology,
and Regents of the University of California. This data license is separate
from the repository's MIT software license. The SenMayo source remains Saul
et al. 2022 as described in the README. The historical Fridman and SASP custom
panels remain explicitly unverified; they are not substituted for these
source-pinned directional sets.

The signed method computes independent control-subtracted UP and DOWN module
scores, then subtracts DOWN from UP without z-standardization. It follows the
database direction but is a project scoring procedure. Neither it nor the
synthetic specificity checks establish biological senescence or rejuvenation.

## Public expression check (2026-10-04)

[GEO GSE268487](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268487)
supplies the actual human LF1 fibroblast count matrix for the first external
check. Three deposited RNA-seq libraries per condition are labelled as
proliferating, quiescent and senescent, with replicate numbers in GSM titles.
Donor and experiment block IDs are absent and are kept null. The source record
provides no dataset-specific license; public accessibility is not an MIT grant.
The raw counts and full metadata are downloaded outside git.

The linked [Dalgarno et al. 2025 preprint](https://pmc.ncbi.nlm.nih.gov/articles/PMC12262444/)
has DOI 10.1101/2025.06.11.659151 and a CC BY-NC-ND 4.0 article license. This
repository links to it and does not redistribute its text or figures. A
post-SenMayo accession supports an external study check, not verified independent
donors/laboratories or independent functional assay validation. LF1 is also used
in older senescence studies; related models do not establish independent donors.

Ensembl identifiers map through the actual HGNC complete-set snapshot retrieved
on 2026-10-04. The committed gzip snapshot contains approved, unambiguous
Ensembl-to-symbol pairs derived by `reduce_hgnc`; it has no expression values.
Three ambiguous identifiers are omitted. HGNC imposes no reuse restrictions
and requests attribution to **HUGO Gene Nomenclature Committee at the University
of Cambridge**, [genenames.org](https://www.genenames.org/about/), RRID:SCR_002827.
Current approved symbols are used without aliases; published panel members are
not altered, so renamed legacy symbols can remain missing from custom panels.

[sources.json](validation/GSE268487/sources.json) pins exact source URLs, byte
counts, SHA-256 values, independence limits, factual metadata snapshot hashes,
mapping hash and the unchanged prespecified plan hash. The actual count gzip
SHA-256 is `ccf26d8b418f9c701c609195e6da0eb742accd5528e24d8811fb9f1fb6dc159a`.
The [evaluation report](validation/GSE268487/REPORT.md) presents negative as well
as positive results and explicitly records a discovered PAR_Y importer correction.

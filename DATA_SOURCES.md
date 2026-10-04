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

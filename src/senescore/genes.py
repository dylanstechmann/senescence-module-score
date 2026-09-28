"""Human SenMayo symbols from Saul et al., Nature Communications 2022.

DOI: 10.1038/s41467-022-32552-1
The gene set is reproduced so a table can be scored offline.
This file is not the paper, and the score below is not their GSEA procedure.
CDKN1A and CDKN2A are intentionally absent: they were not members of this set.
"""

from __future__ import annotations

SENMAYO_HUMAN = (
    "ACVR1B", "ANG", "ANGPT1", "ANGPTL4", "AREG", "AXL", "BEX3", "BMP2", "BMP6", "C3",
    "CCL1", "CCL13", "CCL16", "CCL2", "CCL20", "CCL24", "CCL26", "CCL3", "CCL3L1", "CCL4",
    "CCL5", "CCL7", "CCL8", "CD55", "CD9", "CSF1", "CSF2", "CSF2RB", "CST4", "CTNNB1",
    "CTSB", "CXCL1", "CXCL10", "CXCL12", "CXCL16", "CXCL2", "CXCL3", "CXCL8", "CXCR2", "DKK1",
    "EDN1", "EGF", "EGFR", "EREG", "ESM1", "ETS2", "FAS", "FGF1", "FGF2", "FGF7",
    "GDF15", "GEM", "GMFG", "HGF", "HMGB1", "ICAM1", "ICAM3", "IGF1", "IGFBP1", "IGFBP2",
    "IGFBP3", "IGFBP4", "IGFBP5", "IGFBP6", "IGFBP7", "IL10", "IL13", "IL15", "IL18", "IL1A",
    "IL1B", "IL2", "IL32", "IL6", "IL6ST", "IL7", "INHA", "IQGAP2", "ITGA2", "ITPKA",
    "JUN", "KITLG", "LCP1", "MIF", "MMP1", "MMP10", "MMP12", "MMP13", "MMP14", "MMP2",
    "MMP3", "MMP9", "NAP1L4", "NRG1", "PAPPA", "PECAM1", "PGF", "PIGF", "PLAT", "PLAU",
    "PLAUR", "PTBP1", "PTGER2", "PTGES", "RPS6KA5", "SCAMP4", "SELPLG", "SEMA3F", "SERPINB4", "SERPINE1",
    "SERPINE2", "SPP1", "SPX", "TIMP2", "TNF", "TNFRSF10C", "TNFRSF11B", "TNFRSF1A", "TNFRSF1B", "TUBGCP2",
    "VEGFA", "VEGFC", "VGF", "WNT16", "WNT2",
)

FRIDMAN_SENESCENCE = (
    "ACTA2", "ANGPTL4", "AXL", "BEX3", "BMP2", "CCND1", "CCNE1", "CDKN1A", "CDKN2A",
    "COL1A1", "COL1A2", "COL3A1", "CSF1", "CTGF", "CXCL1", "CXCL8", "CYR61", "DDB2",
    "EGR1", "ETS2", "FAS", "FN1", "GADD45A", "GDF15", "HIC1", "ICAM1", "IGFBP2",
    "IGFBP3", "IGFBP5", "IGFBP7", "IL1A", "IL1B", "IL6", "ING1", "IRF1", "JUN",
    "LMNB1", "MDM2", "MIF", "MMP1", "MMP10", "MMP2", "MMP3", "MYC", "NDN",
    "NFKB1", "NOX4", "PCNA", "PAPPA", "PLAU", "PLAUR", "PML", "PTBP1", "RPS6KA5",
    "SERPINE1", "SERPINE2", "SPARC", "SPP1", "TGFB1", "THBS1", "TIMP1", "TIMP2",
    "TNFRSF10B", "TP53", "VEGFA", "VIM",
)

SASP_COPPE = (
    "AREG", "CCL1", "CCL2", "CCL20", "CCL26", "CCL3", "CCL4", "CCL5", "CCL7", "CCL8",
    "CSF2", "CXCL1", "CXCL2", "CXCL3", "CXCL8", "CXCL10", "CXCL12", "EGF", "FAS",
    "FGF2", "GDF15", "HGF", "ICAM1", "IGFBP2", "IGFBP3", "IGFBP4", "IGFBP6", "IGFBP7",
    "IL1A", "IL1B", "IL6", "IL6ST", "IL7", "IL13", "IL15", "MIF", "MMP1", "MMP3",
    "MMP10", "MMP12", "MMP13", "MMP14", "PAPPA", "PECAM1", "PGF", "PLAT", "PLAU",
    "PLAUR", "SERPINE1", "TIMP1", "TIMP2", "TNF", "TNFRSF1A", "TNFRSF1B", "VEGFA",
)

CITATION = "Saul et al., Nature Communications 13, 4827 (2022). https://doi.org/10.1038/s41467-022-32552-1"
ORTHOGONAL = ("CDKN1A", "CDKN2A")

GENE_SETS: dict[str, dict] = {
    "senmayo": {
        "name": "SenMayo",
        "symbols": SENMAYO_HUMAN,
        "citation": CITATION,
        "description": "Human SenMayo 125-gene senescence signature (Saul et al. 2022)",
        "orthogonal": ORTHOGONAL,
    },
    "fridman": {
        "name": "Fridman Senescence",
        "symbols": FRIDMAN_SENESCENCE,
        "citation": "Fridman et al., Cell Cycle 5(12), 2854-2860 (2006). https://doi.org/10.4161/cc.5.12.2854",
        "description": "Fridman 66-gene senescence and p53/p21 pathway transcriptional signature",
        "orthogonal": (),
    },
    "sasp": {
        "name": "SASP Coppé",
        "symbols": SASP_COPPE,
        "citation": "Coppé et al., PLoS Biology 6(12), e301 (2008). https://doi.org/10.1371/journal.pbio.0060301",
        "description": "Coppé 2008 core senescence-associated secretory phenotype (SASP) signature",
        "orthogonal": ORTHOGONAL,
    },
}


def get_gene_set(name: str) -> dict:
    key = name.strip().lower()
    if key not in GENE_SETS:
        options = ", ".join(sorted(GENE_SETS.keys()))
        raise KeyError(f"unknown gene set '{name}'; available: {options}")
    return GENE_SETS[key]


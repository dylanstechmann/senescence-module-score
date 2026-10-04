"""Human SenMayo symbols from Saul et al., Nature Communications 2022.

DOI: 10.1038/s41467-022-32552-1
The gene set is reproduced so a table can be scored offline.
This file is not the paper, and the score below is not their GSEA procedure.
CDKN1A and CDKN2A are intentionally absent: they were not members of this set.
"""

from __future__ import annotations

import json
import hashlib
from pathlib import Path

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

# This historical project panel is retained for backwards-compatible scoring,
# but it is not the source-aligned MSigDB Fridman set. In particular it mixes
# genes from the source's up/down signatures, so do not interpret its direction.
FRIDMAN_CUSTOM_PANEL = (
    "ACTA2", "ANGPTL4", "AXL", "BEX3", "BMP2", "CCND1", "CCNE1", "CDKN1A", "CDKN2A",
    "COL1A1", "COL1A2", "COL3A1", "CSF1", "CTGF", "CXCL1", "CXCL8", "CYR61", "DDB2",
    "EGR1", "ETS2", "FAS", "FN1", "GADD45A", "GDF15", "HIC1", "ICAM1", "IGFBP2",
    "IGFBP3", "IGFBP5", "IGFBP7", "IL1A", "IL1B", "IL6", "ING1", "IRF1", "JUN",
    "LMNB1", "MDM2", "MIF", "MMP1", "MMP10", "MMP2", "MMP3", "MYC", "NDN",
    "NFKB1", "NOX4", "PCNA", "PAPPA", "PLAU", "PLAUR", "PML", "PTBP1", "RPS6KA5",
    "SERPINE1", "SERPINE2", "SPARC", "SPP1", "TGFB1", "THBS1", "TIMP1", "TIMP2",
    "TNFRSF10B", "TP53", "VEGFA", "VIM",
)


def _load_msigdbr_asset(filename: str) -> dict:
    """Load a source-pinned MSigDB fixture shipped with the package."""
    path = Path(__file__).resolve().parent / "data" / filename
    entry = json.loads(path.read_text(encoding="utf-8"))
    symbols, entrez_ids = entry.get("symbols"), entry.get("entrez_ids")
    if (not isinstance(symbols, list) or not symbols or len(symbols) != len(set(symbols))
            or not isinstance(entrez_ids, list) or len(entrez_ids) != len(symbols)
            or len(entrez_ids) != len(set(entrez_ids))):
        raise RuntimeError(f"invalid source-pinned MSigDB gene-set asset: {path.name}")
    normalized = json.dumps(
        {"symbols": symbols, "entrez_ids": entrez_ids},
        sort_keys=True, separators=(",", ":"), ensure_ascii=True,
    ).encode("utf-8")
    if entry.get("normalized_membership_sha256") != hashlib.sha256(normalized).hexdigest():
        raise RuntimeError(f"invalid normalized MSigDB membership hash: {path.name}")
    expected = {
        "fridman-senescence-up.json": ("FRIDMAN_SENESCENCE_UP", "M9143", "up", 77),
        "fridman-senescence-down.json": ("FRIDMAN_SENESCENCE_DN", "M9487", "down", 13),
    }[filename]
    if ((entry.get("set_id"), entry.get("systematic_id"), entry.get("direction"), len(symbols))
            != expected or entry.get("species") != "Homo sapiens"):
        raise RuntimeError(f"invalid directional MSigDB metadata: {path.name}")
    return entry


def _validate_fridman_gmt_snapshot() -> None:
    """Verify the directional JSON records are exact members of the pinned GMT."""
    path = Path(__file__).resolve().parent / "data" / "msigdb-c2-v2025.1.Hs-fridman-senescence.gmt"
    lines = path.read_bytes().splitlines(keepends=True)
    expected = {"FRIDMAN_SENESCENCE_UP": _FRIDMAN_UP, "FRIDMAN_SENESCENCE_DN": _FRIDMAN_DOWN}
    if len(lines) != len(expected):
        raise RuntimeError("pinned Fridman GMT snapshot must contain exactly the UP and DN entries")
    seen = set()
    for raw_line in lines:
        fields = raw_line.decode("utf-8").rstrip("\n").split("\t")
        if len(fields) < 3 or fields[0] not in expected or fields[0] in seen:
            raise RuntimeError("pinned Fridman GMT snapshot has invalid or duplicate gene-set entries")
        entry = expected[fields[0]]
        symbols = fields[2:]
        if (symbols != entry["symbols"]
                or hashlib.sha256(raw_line).hexdigest() != entry["source"]["raw_line_sha256"]):
            raise RuntimeError(f"pinned MSigDB GMT line does not match the {fields[0]} JSON record")
        seen.add(fields[0])
    if seen != set(expected):
        raise RuntimeError("pinned Fridman GMT snapshot is missing a directional signature")


_FRIDMAN_UP = _load_msigdbr_asset("fridman-senescence-up.json")
_FRIDMAN_DOWN = _load_msigdbr_asset("fridman-senescence-down.json")
_validate_fridman_gmt_snapshot()
FRIDMAN_SENESCENCE_UP = tuple(_FRIDMAN_UP["symbols"])
FRIDMAN_SENESCENCE_DOWN = tuple(_FRIDMAN_DOWN["symbols"])
FRIDMAN_SIGNED_CITATION = (
    "Fridman AL, Tainsky MA. Oncogene. 2008;27:5975-5987. "
    "DOI:10.1038/onc.2008.213; MSigDB C2 release 2025.1.Hs, Table 2S."
)

# This convenience panel is not a machine-readable transcription of a
# particular supplementary table; describe it as custom until that is sourced.
SASP_CUSTOM_PANEL = (
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
        "source_status": "published_gene_set",
    },
    "fridman_up": {
        "name": "MSigDB FRIDMAN_SENESCENCE_UP",
        "symbols": FRIDMAN_SENESCENCE_UP,
        "citation": FRIDMAN_SIGNED_CITATION,
        "description": "77 source-transcribed genes upregulated in senescent cells (Table 2S)",
        "orthogonal": (),
        "source_status": "published_gene_set",
        "direction": "up",
        "systematic_id": _FRIDMAN_UP["systematic_id"],
        "source": _FRIDMAN_UP["source"],
        "normalized_membership_sha256": _FRIDMAN_UP["normalized_membership_sha256"],
    },
    "fridman_down": {
        "name": "MSigDB FRIDMAN_SENESCENCE_DN",
        "symbols": FRIDMAN_SENESCENCE_DOWN,
        "citation": FRIDMAN_SIGNED_CITATION,
        "description": "13 source-transcribed genes downregulated in senescent cells (Table 2S)",
        "orthogonal": (),
        "source_status": "published_gene_set",
        "direction": "down",
        "systematic_id": _FRIDMAN_DOWN["systematic_id"],
        "source": _FRIDMAN_DOWN["source"],
        "normalized_membership_sha256": _FRIDMAN_DOWN["normalized_membership_sha256"],
    },
    "fridman": {
        "name": "Fridman Custom Panel",
        "symbols": FRIDMAN_CUSTOM_PANEL,
        "citation": "Project-curated custom 66-gene panel; membership and direction not verified against a canonical Fridman signature.",
        "description": "Custom panel inspired by Fridman senescence literature; not the canonical MSigDB FRIDMAN_SENESCENCE_UP set",
        "orthogonal": (),
        "source_status": "custom_unverified",
    },
    "sasp": {
        "name": "SASP Custom Panel",
        "symbols": SASP_CUSTOM_PANEL,
        "citation": "Project-curated custom SASP panel; membership not verified against a specific Coppé supplementary table.",
        "description": "Custom SASP-oriented panel; not a source-transcribed Coppé signature",
        "orthogonal": ORTHOGONAL,
        "source_status": "custom_unverified",
    },
}

# Backwards-compatible names for callers that imported the original constants.
FRIDMAN_SENESCENCE = FRIDMAN_CUSTOM_PANEL
SASP_COPPE = SASP_CUSTOM_PANEL


def get_gene_set(name: str) -> dict:
    key = name.strip().lower()
    if key not in GENE_SETS:
        options = ", ".join(sorted(GENE_SETS.keys()))
        raise KeyError(f"unknown gene set '{name}'; available: {options}")
    return GENE_SETS[key]

#!/usr/bin/env python3
"""
Build a reviewer CSV of LOW-confidence CLIP samples with GSM-backed naming errors.

Flags likely sample-name inaccuracies (cell line, protein target, replicate/section
numbering) using GEO GSM title, characteristics, and library name as evidence.

Example:
  python3 flow_sample_name_review_flags.py \\
    --audit-csv projects/flow_sample_name_audit_all.csv \\
    --gsm-cache projects/flow_gsm_cache.json \\
    --output-csv projects/flow_sample_name_review_flags.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from flow_sample_name_audit import PROTEIN_ALIASES, GeoGsm, cell_matches, norm_token, protein_matches
from flow_sample_name_audit_all import library_name_from_geo, parse_geo_catalog, split_protein_catalog

FLOW_SAMPLE_URL = "https://app.flow.bio/samples/{sample_id}"
GSM_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}"

CELL_LINES = {"A431", "HEK293T", "HEK293", "K562", "HepG2", "HCT116", "mESC", "nESC", "N2A", "ESC"}
CELL_EQUIV = {
    "nesc": {"nesc", "esc", "mesc"},
    "esc": {"nesc", "esc", "mesc"},
    "hek293": {"hek293", "hek293t"},
    "hek293t": {"hek293", "hek293t"},
    "n2a": {"n2a", "neuro2a"},
}
CLONE_TOKEN_RE = re.compile(r"^c\d+$", re.I)

# Common gene symbol aliases not always in the base audit table
EXTRA_PROTEIN_ALIASES = {
    "TARDBP": ["TARDBP", "TDP43"],
    "GARS1": ["GARS1", "GARS"],
    "FMR1": ["FMR1", "FMRP"],
}


def protein_matches_extended(observed: str, expected: str) -> Tuple[float, str]:
    score, note = protein_matches(observed, expected)
    if score >= 1.0:
        return score, note
    o, e = norm_token(observed), norm_token(expected)
    for _canon, aliases in EXTRA_PROTEIN_ALIASES.items():
        alias_norm = {norm_token(a) for a in aliases}
        if o in alias_norm and e in alias_norm:
            return 1.0, "alias"
    return score, note


@dataclass
class NamingFlag:
    error_type: str
    name_observed: str
    gsm_expected: str
    evidence: str
    severity: str  # high | medium


def load_gsm_cache(path: Path) -> Dict[str, GeoGsm]:
    if not path.is_file():
        return {}
    data = json.loads(path.read_text(encoding="utf-8"))
    return {k: GeoGsm(**v) for k, v in data.items()}


def extract_geo_protein_pair(title: str) -> Tuple[str, str]:
    parts = [p.strip() for p in title.split(",") if p.strip()]
    primary = parts[0] if parts else ""
    secondary = ""
    for p in parts[1:]:
        pl = p.lower()
        if pl in {"irclipv2", "re-clip", "iclip", "iclip2"} or "egf" in pl or pl.startswith("rep"):
            break
        if p in CELL_LINES or re.fullmatch(r"\d+min EGF", p, re.I):
            break
        secondary = p
        break
    return primary, secondary


def extract_geo_protein(title: str, row: Dict[str, str], library_name: str = "") -> str:
    if not title:
        return ""
    patterns = [
        r"^([^(,]+?)\s*\(\d+\)\s+eCLIP",
        r"^([A-Za-z0-9]+),\s*irCLIPv2",
        r"^([A-Za-z0-9]+),\s*Re-CLIP",
        r"^([A-Za-z0-9]+),\s*iCLIP",
        r"iCLIP-seq of (?:Flag-)?([A-Za-z0-9]+)",
    ]
    for pat in patterns:
        m = re.search(pat, title, re.I)
        if m:
            return m.group(1).strip()
    if library_name:
        lib_tok = library_name.split("_")[0]
        if lib_tok and not lib_tok.lower().startswith("control"):
            return lib_tok
    flow_target = (row.get("flow_purification_target") or "").strip()
    if flow_target and norm_token(flow_target) not in {"sminput", "input", "igg"}:
        return flow_target
    m = re.search(r"^([A-Za-z0-9]+),", title)
    if m:
        return m.group(1).strip()
    return ""


def extract_geo_cell_title(title: str) -> str:
    for pat in [r"in\s+(\w+)\s+cells", r"from\s+(\w+)\s*\(", r"from\s+(\w+)\s+cells"]:
        m = re.search(pat, title, re.I)
        if m:
            return m.group(1)
    parts = [p.strip() for p in title.split(",") if p.strip()]
    for p in parts:
        if p in CELL_LINES:
            return p
    return ""


def extract_geo_rep(title: str, library_name: str = "") -> str:
    for src in (title, library_name):
        if not src:
            continue
        m = re.search(r"(?:replicate|rep)\s*(\d+)", src, re.I)
        if m:
            return m.group(1)
        m = re.search(r"\bR(\d+)\b", src)
        if m:
            return m.group(1)
    return ""


def extract_geo_section(title: str) -> str:
    m = re.search(r"Section\s*(\d+)", title, re.I)
    return m.group(1) if m else ""


def parse_irclip_name(name: str) -> Dict[str, str]:
    tail = re.match(
        r"^(?P<head>.+?)(?:_(?P<rep>R\d+))?(?:_(?P<section>S\d+))?_(?P<date>\d{8})$",
        name,
    )
    if not tail:
        return {}
    head = tail.group("head")
    core = re.match(
        r"^(?P<prefix>.+)_Hs_(?P<cell>[^_]+)_(?P<method>[^_]+)_(?P<condition>.+)$",
        head,
    )
    if not core:
        return {}
    prefix_tokens = core.group("prefix").split("_")
    return {
        "format": "irclip",
        "protein_primary": prefix_tokens[0],
        "protein_secondary": prefix_tokens[1] if len(prefix_tokens) > 1 else "",
        "cell": core.group("cell"),
        "rep": (tail.group("rep") or "").lstrip("R"),
        "section": (tail.group("section") or "").lstrip("S"),
    }


def parse_encore_name(name: str) -> Dict[str, str]:
    m = re.match(
        r"^(?P<protein_token>[A-Za-z0-9]+)_Hs_(?P<cell>[^_]+)_(?P<condition>INPUT|IP)_rep(?P<rep>\d+)$",
        name,
        re.I,
    )
    if not m:
        return {}
    return {
        "format": "encore",
        "protein_token": m.group("protein_token"),
        "cell": m.group("cell"),
        "condition": m.group("condition").upper(),
        "rep": m.group("rep"),
    }


def cells_equivalent(observed: str, expected: str) -> bool:
    o, e = norm_token(observed), norm_token(expected)
    if not o or not e:
        return False
    if o == e:
        return True
    for variants in CELL_EQUIV.values():
        if o in variants and e in variants:
            return True
    return cell_matches(observed, expected)[0] >= 1.0


KNOWN_CELL_TOKENS = {
    "a431", "hek293", "hek293t", "k562", "hepg2", "hct116", "mesc", "nesc", "n2a", "esc",
    "mdamb231", "mdamb436", "mcf10a", "sum149", "cvb", "h1motor", "hepg2",
}


def looks_like_cell_token(token: str) -> bool:
    if not token or token.lower() in {"input", "ip", "flag", "control", "rep1", "rep2"}:
        return False
    if CLONE_TOKEN_RE.fullmatch(token) or re.fullmatch(r"p\d+h", token, re.I):
        return False
    return norm_token(token) in KNOWN_CELL_TOKENS or norm_token(token) in {norm_token(c) for c in CELL_LINES}


def extract_cell_token(name: str) -> str:
    m = re.search(r"_(?P<cell>[^_]+)_Mm_", name)
    if m and looks_like_cell_token(m.group("cell")):
        return m.group("cell")
    m = re.search(r"_Mm_(?P<cell>[^_]+)", name)
    if m:
        tok = m.group("cell")
        if tok.lower() not in {"input", "ip", "ctfusion"} and not re.fullmatch(r"rep\d+", tok, re.I):
            return tok
    m = re.search(r"_(?P<cell>[^_]+)_Hs_", name)
    if m and looks_like_cell_token(m.group("cell")):
        return m.group("cell")
    m = re.search(r"_Hs_(?P<cell>[^_]+)(?:_|$)", name)
    if m and looks_like_cell_token(m.group("cell")):
        return m.group("cell")
    if m:
        return m.group("cell")
    m = re.search(r"_(?P<cell>[A-Za-z0-9][A-Za-z0-9-]*)_Hs_rep(?P<rep>\d+)", name, re.I)
    if m:
        return m.group("cell")
    return ""


def parse_generic_name(name: str) -> Dict[str, str]:
    tokens = name.split("_") if name else []
    out: Dict[str, str] = {"format": "generic", "protein_primary": tokens[0] if tokens else ""}
    if len(tokens) >= 2 and tokens[0].upper() == "DB21":
        out["protein_primary"] = tokens[1]
    cell = extract_cell_token(name)
    if cell:
        out["cell"] = cell
    m = re.search(r"_(?P<cell>[A-Za-z0-9][A-Za-z0-9-]*)_Hs_rep(?P<rep>\d+)", name, re.I)
    if m:
        out["cell"] = m.group("cell")
        out["rep"] = m.group("rep")
    if "rep" not in out:
        m = re.search(r"(?:^|_)(?:rep|R)(?P<rep>\d+)(?:_|$)", name, re.I)
        if m:
            out["rep"] = m.group(1)
    m = re.search(r"(?:^|_)S(?P<section>\d+)(?:_|$)", name)
    if m:
        out["section"] = m.group("section")
    return out


def parse_name_tokens(name: str) -> Dict[str, str]:
    for parser in (parse_irclip_name, parse_encore_name):
        parsed = parser(name)
        if parsed:
            return parsed
    return parse_generic_name(name)


def gsm_authoritative_cell(geo: GeoGsm) -> str:
    for candidate in (geo.cell_line, geo.source_name, extract_geo_cell_title(geo.title)):
        if candidate and looks_like_cell_token(candidate.split()[0]):
            return candidate.split()[0]
    m = re.match(r"^(HEK293T?|K562|A431|HepG2|HCT116|HEK293)", geo.title, re.I)
    if m:
        return m.group(1)
    return geo.cell_line or geo.source_name or extract_geo_cell_title(geo.title)


def flow_metadata_agrees_with_gsm(row: Dict[str, str], geo: GeoGsm, field: str, gsm_value: str) -> bool:
    if not gsm_value:
        return False
    if field == "cell":
        flow_val = row.get("flow_source", "")
        return cell_matches(flow_val, gsm_value)[0] >= 1.0
    if field == "protein":
        flow_val = row.get("flow_purification_target", "")
        if norm_token(row.get("flow_condition", "")) == "input" or norm_token(flow_val) == "sminput":
            return True
        return protein_matches_extended(flow_val, gsm_value)[0] >= 1.0
    return False


def detect_naming_flags(row: Dict[str, str], geo: GeoGsm) -> List[NamingFlag]:
    flags: List[NamingFlag] = []
    name = row.get("sample_name", "")
    parsed = parse_name_tokens(name)
    if not parsed:
        return flags

    lib = library_name_from_geo(geo)
    catalog = parse_geo_catalog(geo.title)
    gsm_cell = gsm_authoritative_cell(geo)
    gsm_prot = extract_geo_protein(geo.title, row, lib)
    gsm_rep = extract_geo_rep(geo.title, lib)
    gsm_section = extract_geo_section(geo.title)
    title_cell = extract_geo_cell_title(geo.title)

    # --- cell line in sample name ---
    name_cell = parsed.get("cell", "")
    if name_cell and (
        name_cell.lower() in {"input", "ip", "rep", "control"}
        or CLONE_TOKEN_RE.fullmatch(name_cell)
        or re.fullmatch(r"p\d+h", name_cell, re.I)
    ):
        name_cell = ""
    if name_cell and gsm_cell:
        if cells_equivalent(name_cell, gsm_cell):
            name_cell = ""
    if name_cell and gsm_cell:
        score, _ = cell_matches(name_cell, gsm_cell)
        if score < 1.0 and looks_like_cell_token(name_cell):
            flow_agrees = flow_metadata_agrees_with_gsm(row, geo, "cell", gsm_cell)
            title_matches_name = title_cell and cell_matches(name_cell, title_cell)[0] >= 1.0
            char_conflicts_title = (
                title_cell
                and geo.cell_line
                and cell_matches(title_cell, geo.cell_line)[0] < 1.0
            )
            if flow_agrees or (title_matches_name and char_conflicts_title):
                evidence_parts = [f'GSM title: "{geo.title}"']
                if geo.cell_line:
                    evidence_parts.append(f'GSM cell line characteristic: "{geo.cell_line}"')
                if geo.source_name:
                    evidence_parts.append(f'GSM source name: "{geo.source_name}"')
                if lib:
                    evidence_parts.append(f'GSM library name: "{lib}"')
                if flow_agrees:
                    evidence_parts.append(f'Flow source agrees with GSM: "{row.get("flow_source", "")}"')
                if title_matches_name and char_conflicts_title:
                    evidence_parts.append(
                        f'Name cell "{name_cell}" matches GSM title but characteristic says "{gsm_cell}"'
                    )
                flags.append(
                    NamingFlag(
                        error_type="cell_line",
                        name_observed=name_cell,
                        gsm_expected=gsm_cell,
                        evidence=" | ".join(evidence_parts),
                        severity="high",
                    )
                )

    # --- protein target in sample name ---
    gsm_primary, gsm_secondary = extract_geo_protein_pair(geo.title)
    if parsed.get("format") == "encore":
        name_prot, name_cat = split_protein_catalog(parsed.get("protein_token", ""), catalog)
    elif parsed.get("format") == "irclip" and parsed.get("protein_secondary"):
        # Re-CLIP names encode scaffold_IPtarget; compare IP target token to GSM/Flow.
        name_prot = parsed.get("protein_secondary", "")
        gsm_prot = gsm_secondary or gsm_prot or gsm_primary
    else:
        name_prot = parsed.get("protein_primary", "")
        name_cat = ""

    if name_prot and gsm_prot and norm_token(name_prot) not in {"control", "input", "igg"}:
        p_score, _ = protein_matches_extended(name_prot, gsm_prot)
        flow_target = norm_token(row.get("flow_purification_target", ""))
        name_contains_target = (
            norm_token(gsm_prot) in norm_token(name_prot)
            or (flow_target and flow_target in norm_token(name_prot))
            or (flow_target and norm_token(name_prot).startswith(flow_target))
        )
        if p_score < 1.0 and name_contains_target:
            p_score = 1.0
        if p_score < 1.0 and flow_metadata_agrees_with_gsm(row, geo, "protein", gsm_prot):
            evidence = f'GSM title: "{geo.title}"'
            if lib:
                evidence += f' | GSM library name: "{lib}"'
            evidence += f' | Flow purification target agrees: "{row.get("flow_purification_target", "")}"'
            flags.append(
                NamingFlag(
                    error_type="protein_target",
                    name_observed=name_prot,
                    gsm_expected=gsm_prot,
                    evidence=evidence,
                    severity="high",
                )
            )

    # --- replicate numbering ---
    name_rep = parsed.get("rep", "")
    if name_rep and gsm_rep and name_rep != gsm_rep:
        evidence = f'GSM title: "{geo.title}"'
        if lib:
            evidence += f' | GSM library name: "{lib}"'
        flags.append(
            NamingFlag(
                error_type="replicate_number",
                name_observed=name_rep,
                gsm_expected=gsm_rep,
                evidence=evidence,
                severity="high",
            )
        )

    # --- section numbering (irCLIP molecular-weight fractions) ---
    name_section = parsed.get("section", "")
    if name_section and gsm_section and name_section != gsm_section:
        flags.append(
            NamingFlag(
                error_type="section_number",
                name_observed=name_section,
                gsm_expected=gsm_section,
                evidence=f'GSM title: "{geo.title}"',
                severity="high",
            )
        )

    # --- INPUT/IP condition token ---
    if parsed.get("format") == "encore":
        geo_cond = "INPUT" if re.search(r"\bINPUT\b", geo.title, re.I) else ""
        if re.search(r"\bIP\b", geo.title, re.I):
            geo_cond = "IP"
        name_cond = parsed.get("condition", "")
        if geo_cond and name_cond and norm_token(name_cond) != norm_token(geo_cond):
            flags.append(
                NamingFlag(
                    error_type="condition_token",
                    name_observed=name_cond,
                    gsm_expected=geo_cond,
                    evidence=f'GSM title: "{geo.title}"',
                    severity="high",
                )
            )

    return flags


def _suggested_issue(flags: List[NamingFlag]) -> str:
    parts: List[str] = []
    for f in flags:
        if f.error_type == "cell_line":
            parts.append(f'Sample name cell token "{f.name_observed}" may be wrong; GSM expects "{f.gsm_expected}"')
        elif f.error_type == "protein_target":
            parts.append(f'Sample name protein token "{f.name_observed}" may be wrong; GSM/Flow expect "{f.gsm_expected}"')
        elif f.error_type == "replicate_number":
            parts.append(f'Sample name replicate "{f.name_observed}" may be wrong; GSM expects rep "{f.gsm_expected}"')
        elif f.error_type == "section_number":
            parts.append(f'Sample name section "{f.name_observed}" may be wrong; GSM expects section "{f.gsm_expected}"')
        elif f.error_type == "condition_token":
            parts.append(f'Sample name condition "{f.name_observed}" may be wrong; GSM expects "{f.gsm_expected}"')
    return " | ".join(parts)


def flag_row_to_dict(row: Dict[str, str], geo: GeoGsm, flags: List[NamingFlag]) -> Dict[str, Any]:
    primary = flags[0]
    types = sorted({f.error_type for f in flags})
    severities = sorted({f.severity for f in flags})
    high_flags = [f for f in flags if f.severity == "high"]

    return {
        "sample_id": row.get("sample_id", ""),
        "sample_name": row.get("sample_name", ""),
        "project_id": row.get("project_id", ""),
        "project_name": row.get("project_name", ""),
        "flow_sample_url": FLOW_SAMPLE_URL.format(sample_id=row.get("sample_id", "")),
        "gsm": row.get("lookup_accession", ""),
        "gsm_url": GSM_URL.format(gsm=row.get("lookup_accession", "")),
        "confidence_score": row.get("confidence_score", ""),
        "confidence_band": row.get("confidence_band", ""),
        "naming_error_types": ";".join(types),
        "flag_severity": "high" if "high" in severities else "medium",
        "flag_count": len(flags),
        "primary_error_type": primary.error_type,
        "name_value_observed": primary.name_observed,
        "gsm_value_expected": primary.gsm_expected,
        "gsm_evidence": primary.evidence,
        "all_naming_flags": " || ".join(
            f"{f.error_type}:{f.name_observed}->{f.gsm_expected} [{f.severity}] ({f.evidence})" for f in flags
        ),
        "gsm_title": geo.title,
        "gsm_cell_line_characteristic": geo.cell_line,
        "gsm_source_name": geo.source_name,
        "gsm_library_name": library_name_from_geo(geo),
        "flow_source": row.get("flow_source", ""),
        "flow_purification_target": row.get("flow_purification_target", ""),
        "flow_condition": row.get("flow_condition", ""),
        "flow_method": row.get("flow_method", ""),
        "metadata_agrees_with_gsm": "yes" if high_flags else "no",
        "review_priority": "actionable" if high_flags else "review",
        "suggested_issue": _suggested_issue(flags),
    }


def flags_from_audit_field_scores(row: Dict[str, str], geo: GeoGsm) -> List[NamingFlag]:
    """Supplement parser-based flags using per-field audit scores."""
    flags: List[NamingFlag] = []
    try:
        field_scores = json.loads(row.get("field_scores_json") or "[]")
    except json.JSONDecodeError:
        return flags
    by_field = {fs["field"]: fs for fs in field_scores}
    lib = library_name_from_geo(geo)
    parsed = parse_name_tokens(row.get("sample_name", ""))

    cell_fs = by_field.get("cell")
    if cell_fs and cell_fs.get("score") == 0.0:
        name_cell = parsed.get("cell", "")
        gsm_cell = gsm_authoritative_cell(geo)
        if name_cell and gsm_cell and not cells_equivalent(name_cell, gsm_cell) and looks_like_cell_token(name_cell):
            flow_agrees = flow_metadata_agrees_with_gsm(row, geo, "cell", gsm_cell)
            flags.append(
                NamingFlag(
                    error_type="cell_line",
                    name_observed=name_cell,
                    gsm_expected=gsm_cell,
                    evidence=(
                        f'GSM title: "{geo.title}"'
                        + (f' | GSM cell line: "{geo.cell_line}"' if geo.cell_line else "")
                        + f' | Flow source: "{row.get("flow_source", "")}"'
                        + (" | Flow source matches GSM" if flow_agrees else "")
                    ),
                    severity="high" if flow_agrees else "medium",
                )
            )

    prot_fs = by_field.get("protein")
    if prot_fs and prot_fs.get("score") == 0.0:
        gsm_prot = extract_geo_protein(geo.title, row, lib)
        name_prot = parsed.get("protein_primary", "")
        if parsed.get("format") == "irclip" and parsed.get("protein_secondary"):
            name_prot = parsed.get("protein_secondary", name_prot)
        flow_target = norm_token(row.get("flow_purification_target", ""))
        if name_prot and gsm_prot and protein_matches_extended(name_prot, gsm_prot)[0] < 1.0:
            if not (flow_target and (flow_target in norm_token(name_prot) or norm_token(name_prot).startswith(flow_target))):
                flow_agrees = flow_metadata_agrees_with_gsm(row, geo, "protein", gsm_prot)
                flags.append(
                    NamingFlag(
                        error_type="protein_target",
                        name_observed=name_prot,
                        gsm_expected=gsm_prot,
                        evidence=(
                            f'GSM title: "{geo.title}"'
                            + (f' | GSM library name: "{lib}"' if lib else "")
                            + f' | Flow target: "{row.get("flow_purification_target", "")}"'
                        ),
                        severity="high" if flow_agrees else "medium",
                    )
                )
    return flags


def build_review_rows(
    audit_rows: Sequence[Dict[str, str]],
    gsm_cache: Dict[str, GeoGsm],
    *,
    low_only: bool = True,
    gsm_only: bool = True,
    include_actionable_non_low: bool = True,
) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for row in audit_rows:
        if gsm_only and row.get("lookup_mode") != "gsm":
            continue
        gsm = row.get("lookup_accession", "")
        geo = gsm_cache.get(gsm)
        if not geo or geo.fetch_error:
            continue
        flags = detect_naming_flags(row, geo)
        extra = flags_from_audit_field_scores(row, geo)
        seen = {(f.error_type, f.name_observed, f.gsm_expected) for f in flags}
        for f in extra:
            key = (f.error_type, f.name_observed, f.gsm_expected)
            if key not in seen:
                flags.append(f)
                seen.add(key)
        if not flags:
            continue
        review = flag_row_to_dict(row, geo, flags)
        is_low = row.get("confidence_band") == "LOW"
        is_actionable = review["review_priority"] == "actionable"
        if low_only and not is_low and not (include_actionable_non_low and is_actionable):
            continue
        out.append(review)

    out.sort(
        key=lambda r: (
            0 if r["flag_severity"] == "high" else 1,
            -int(r["flag_count"]),
            float(r["confidence_score"] or 0),
            r["project_name"],
            r["sample_name"],
        )
    )
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--audit-csv", default="/home/mikej10/advbfx/projects/flow_sample_name_audit_all.csv")
    ap.add_argument("--gsm-cache", default="/home/mikej10/advbfx/projects/flow_gsm_cache.json")
    ap.add_argument("--output-csv", default="/home/mikej10/advbfx/projects/flow_sample_name_review_flags.csv")
    ap.add_argument("--include-medium", action="store_true", help="Also include MEDIUM band samples with naming flags")
    ap.add_argument(
        "--actionable-only",
        action="store_true",
        help="Only high-severity flags where Flow metadata corroborates GSM",
    )
    ap.add_argument(
        "--all-flags",
        action="store_true",
        help="Include all medium-severity naming flags (no extra filtering)",
    )
    args = ap.parse_args()
    actionable_only = args.actionable_only
    broad_review = not args.actionable_only and not args.all_flags

    with open(args.audit_csv, newline="", encoding="utf-8") as f:
        audit_rows = list(csv.DictReader(f))
    gsm_cache = load_gsm_cache(Path(args.gsm_cache))

    review_rows = build_review_rows(
        audit_rows,
        gsm_cache,
        low_only=not args.include_medium,
        gsm_only=True,
    )
    if actionable_only:
        review_rows = [r for r in review_rows if r["flag_severity"] == "high"]
    elif broad_review:
        review_rows = [
            r
            for r in review_rows
            if r["flag_severity"] == "high"
            or (
                r["flag_severity"] == "medium"
                and r["confidence_band"] == "LOW"
                and any(
                    t in r["naming_error_types"].split(";")
                    for t in ("cell_line", "protein_target", "replicate_number", "section_number")
                )
            )
        ]

    if review_rows:
        with open(args.output_csv, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(review_rows[0].keys()))
            w.writeheader()
            w.writerows(review_rows)

    summary = {
        "audit_rows": len(audit_rows),
        "flagged_rows": len(review_rows),
        "high_severity": sum(1 for r in review_rows if r["flag_severity"] == "high"),
        "by_error_type": {},
        "by_project": {},
    }
    for r in review_rows:
        for et in r["naming_error_types"].split(";"):
            if et:
                summary["by_error_type"][et] = summary["by_error_type"].get(et, 0) + 1
        pn = r["project_name"]
        summary["by_project"][pn] = summary["by_project"].get(pn, 0) + 1

    summary_path = Path(args.output_csv).with_suffix(".summary.json")
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(f"Wrote {args.output_csv} ({len(review_rows)} flagged samples)")
    print(f"Summary: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

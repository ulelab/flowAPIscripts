#!/usr/bin/env python3
"""
Audit Flow.bio sample names against GEO GSM metadata.

Fetches GSM records (GEO acc.cgi text view), parses comma-separated titles and
characteristics, compares each underscore-delimited name token to expected values,
and emits per-field scores plus an overall confidence score (0–100).

Example:
  python3 flow_sample_name_audit.py \\
    --input-csv /path/to/project_samples.csv \\
    --output-csv /path/to/sample_name_audit.csv \\
    --summary-json /path/to/sample_name_audit_summary.json
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import re
import time
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

# Canonical protein symbols and GEO title aliases (HuR, THOC4, etc.)
PROTEIN_ALIASES: Dict[str, List[str]] = {
    "ALYREF": ["ALYREF", "THOC4"],
    "ELAVL1": ["ELAVL1", "HUR"],
    "FUS": ["FUS"],
    "HNRNPC": ["HNRNPC"],
    "HNRNPM": ["HNRNPM"],
    "HNRNPU": ["HNRNPU"],
    "HNRNPFH": ["HNRNPFH", "HNRNPHF", "HNRNPHF", "HNRNPFH"],
    "HNRNPHF": ["HNRNPHF", "HNRNPFH", "HNRNPHF", "HNRNPFH"],
    "IGG": ["IGG"],
    "KHSRP": ["KHSRP"],
    "UPF1": ["UPF1"],
}

CELL_LINES = {"A431", "HEK293T"}
SPECIES_CODES = {"HS": "Hs", "MM": "Mm", "GG": "Gg"}

METHOD_ALIASES = {
    "irclip": ["irclip", "irclipv2"],
    "reclip": ["reclip", "re-clip"],
}

FIELD_WEIGHTS = {
    "protein": 25.0,
    "cell": 25.0,
    "condition": 20.0,
    "method": 15.0,
    "replicate": 8.0,
    "section": 4.0,
    "species": 3.0,
}


def norm_token(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", (s or "").strip().lower())


def protein_matches(observed: str, expected: str) -> Tuple[float, str]:
    o = norm_token(observed)
    e = norm_token(expected)
    if not o or not e:
        return 0.0, "missing"
    if o == e:
        return 1.0, "exact"
    for canon, aliases in PROTEIN_ALIASES.items():
        alias_norm = {norm_token(a) for a in aliases}
        if o in alias_norm and e in alias_norm:
            return 1.0, f"alias:{canon}"
    return 0.0, "mismatch"


def cell_matches(observed: str, expected: str) -> Tuple[float, str]:
    o = norm_token(observed)
    e = norm_token(expected)
    if not o or not e:
        return 0.0, "missing"
    if o == e:
        return 1.0, "exact"
    if o == "hek293" and e == "hek293t":
        return 1.0, "hek293=hek293t"
    if o == "hek293t" and e == "hek293":
        return 1.0, "hek293t=hek293"
    return 0.0, "mismatch"


def method_matches(observed: str, expected: str) -> Tuple[float, str]:
    o = norm_token(observed)
    e = norm_token(expected)
    if not o or not e:
        return 0.0, "missing"
    for _key, aliases in METHOD_ALIASES.items():
        if o in aliases and e in aliases:
            return 1.0, "exact"
    return 0.0, "mismatch"


def condition_from_name_token(token: str) -> str:
    t = token.strip()
    mapping = {
        "EGF0": "0min EGF",
        "EGF15": "15min EGF",
        "EGF30": "30min EGF",
        "EGF60": "60min EGF",
        "UVC": "UVC",
        "NOUV": "noUV",
        "110-350UV": "110-350UV",
        "140-350UV": "140-350UV",
        "ALLUV": "allUV",
    }
    return mapping.get(t.upper(), t)


def condition_matches(name_token: str, flow_condition: str, geo_title: str, geo_treatment: str) -> Tuple[float, str]:
    name_cond = condition_from_name_token(name_token)
    flow_cond = (flow_condition or "").strip()
    title_low = geo_title.lower()
    treat_low = geo_treatment.lower()

    candidates: List[str] = []
    if flow_cond:
        candidates.append(flow_cond)
    if "nouv" in title_low or "no uv" in title_low:
        candidates.append("noUV")
    if "0min egf" in title_low:
        candidates.append("0min EGF")
    if "15min egf" in title_low:
        candidates.append("15min EGF")
    if "30min egf" in title_low:
        candidates.append("30min EGF")
    if "60min egf" in title_low:
        candidates.append("60min EGF")
    if "110-350uv" in title_low:
        candidates.append("110-350UV")
    if "140-350uv" in title_low:
        candidates.append("140-350UV")
    if "uv cross-linked" in treat_low or "uv crosslinked" in treat_low:
        if name_cond == "UVC":
            candidates.append("UVC")
        if name_cond in {"110-350UV", "140-350UV", "allUV"}:
            candidates.append(name_cond)

    if name_cond == "UVC" and not flow_cond and "section" in title_low:
        candidates.append("UVC")

    if not candidates:
        return 0.5, "no_geo_condition"

    for c in candidates:
        if norm_token(c) == norm_token(name_cond):
            return 1.0, f"match:{c}"
        # EGF token compact vs expanded
        if norm_token(c).replace("minegf", "") == norm_token(name_cond).replace("egf", ""):
            return 1.0, f"match:{c}"
    return 0.0, f"mismatch:expected_one_of={candidates}"


def modifier_matches(prefix_tokens: List[str], geo_title: str) -> Tuple[float, str]:
    title_low = geo_title.lower()
    modifiers = [t for t in prefix_tokens if t.lower().startswith("si")]
    if not modifiers:
        if "sirna control" in title_low or "sirnacontrol" in title_low:
            return 0.0, "missing_siRNAcontrol_in_name"
        if "siupf1" in title_low or "sihnrnpc" in title_low:
            return 0.0, "missing_si_modifier_in_name"
        return 1.0, "n/a"
    for mod in modifiers:
        m = norm_token(mod)
        if m in norm_token(title_low):
            return 1.0, mod
        if m == "sirnacontrol" and "sirna control" in title_low:
            return 1.0, mod
    return 0.0, f"mismatch:{modifiers}"


@dataclass
class GeoGsm:
    gsm: str
    title: str = ""
    source_name: str = ""
    cell_line: str = ""
    cell_type: str = ""
    treatment_protocol: str = ""
    description: str = ""
    extract_protocol: str = ""
    organism: str = ""
    fetch_error: str = ""

    @property
    def title_parts(self) -> List[str]:
        return [p.strip() for p in self.title.split(",") if p.strip()]


def fetch_gsm(gsm: str, cache: Dict[str, GeoGsm], delay_s: float) -> GeoGsm:
    if gsm in cache:
        return cache[gsm]
    out = GeoGsm(gsm=gsm)
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}&targ=self&form=text&view=full"
    try:
        text = urllib.request.urlopen(url, timeout=90).read().decode("utf-8", errors="replace")
        for line in text.splitlines():
            if line.startswith("!Sample_title = "):
                out.title = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_source_name_ch1 = "):
                out.source_name = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_organism_ch1 = "):
                out.organism = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_treatment_protocol_ch1 = "):
                out.treatment_protocol = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_description = "):
                out.description = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_extract_protocol_ch1 = "):
                out.extract_protocol = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_characteristics_ch1 = "):
                val = line.split("=", 1)[1].strip()
                if ":" in val:
                    key, v = val.split(":", 1)
                    key = key.strip().lower()
                    v = v.strip()
                    if key == "cell line":
                        out.cell_line = v
                    elif key == "cell type":
                        out.cell_type = v
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        out.fetch_error = str(exc)
    cache[gsm] = out
    if delay_s > 0:
        time.sleep(delay_s)
    return out


def parse_geo_title(geo: GeoGsm) -> Dict[str, str]:
    parts = geo.title_parts
    method = ""
    cell_title = ""
    proteins: List[str] = []
    for p in parts:
        pl = p.lower()
        if pl in ("irclipv2", "re-clip"):
            method = "Re-CLIP" if pl == "re-clip" else "irCLIP"
            break
        proteins.append(p)
    for p in parts:
        if p in CELL_LINES:
            cell_title = p
    rep_m = re.search(r"Rep\s*(\d+)", geo.title, re.I)
    sec_m = re.search(r"Section\s*(\d+)", geo.title, re.I)
    return {
        "geo_protein_primary": proteins[0] if proteins else "",
        "geo_protein_secondary": proteins[1] if len(proteins) > 1 else "",
        "geo_method": method,
        "geo_cell_title": cell_title,
        "geo_cell_authoritative": geo.cell_line or geo.source_name,
        "geo_rep": rep_m.group(1) if rep_m else "",
        "geo_section": sec_m.group(1) if sec_m else "",
    }


def parse_sample_name(name: str) -> Dict[str, str]:
    # Fixed tail: optional rep/section tokens then YYYYMMDD date suffix.
    tail = re.match(
        r"^(?P<head>.+?)(?:_(?P<rep>R\d+))?(?:_(?P<section>S\d+))?_(?P<date>\d{8})$",
        name,
    )
    if not tail:
        return {"parse_error": f"unparsed:{name}"}
    head = tail.group("head")
    core = re.match(
        r"^(?P<prefix>.+)_Hs_(?P<cell>[^_]+)_(?P<method>[^_]+)_(?P<condition>.+)$",
        head,
    )
    if not core:
        return {"parse_error": f"unparsed:{name}"}
    prefix_tokens = core.group("prefix").split("_")
    condition_raw = core.group("condition")
    hmw = ""
    if condition_raw.endswith("_HMW"):
        hmw = "HMW"
        condition_raw = condition_raw[: -len("_HMW")]
    return {
        "prefix": core.group("prefix"),
        "prefix_tokens": prefix_tokens,
        "protein_primary": prefix_tokens[0],
        "protein_secondary": prefix_tokens[1] if len(prefix_tokens) > 1 else "",
        "species": "Hs",
        "cell": core.group("cell"),
        "method": core.group("method"),
        "condition_token": condition_raw,
        "hmw": hmw,
        "rep": (tail.group("rep") or "").lstrip("R"),
        "section": (tail.group("section") or "").lstrip("S"),
        "date": tail.group("date"),
    }


@dataclass
class FieldScore:
    field: str
    weight: float
    score: float
    status: str
    observed: str = ""
    expected: str = ""
    note: str = ""


def score_sample(row: Dict[str, str], geo: GeoGsm, parsed: Dict[str, str], geo_parsed: Dict[str, str]) -> Tuple[List[FieldScore], float, str]:
    if parsed.get("parse_error"):
        return [FieldScore("parse", 100.0, 0.0, "fail", note=parsed["parse_error"])], 0.0, "PARSE_ERROR"

    scores: List[FieldScore] = []
    prefix_tokens: List[str] = parsed.get("prefix_tokens", [])

    # Protein primary (IP target for irCLIP; reCLIP scaffold for Re-CLIP)
    if norm_token(parsed["method"]) in METHOD_ALIASES["reclip"]:
        exp_prot = geo_parsed["geo_protein_primary"]
        obs_prot = parsed["protein_primary"]
    else:
        exp_prot = geo_parsed["geo_protein_primary"]
        obs_prot = parsed["protein_primary"]
    p_score, p_note = protein_matches(obs_prot, exp_prot)
    scores.append(
        FieldScore(
            "protein",
            FIELD_WEIGHTS["protein"],
            p_score,
            "pass" if p_score >= 1.0 else "fail",
            observed=obs_prot,
            expected=exp_prot,
            note=p_note,
        )
    )

    # Secondary protein (Re-CLIP second target, siRNA modifiers, IgG)
    if parsed["protein_secondary"]:
        exp_sec = geo_parsed["geo_protein_secondary"]
        obs_sec = parsed["protein_secondary"]
        if norm_token(obs_sec).startswith("si"):
            s_score, s_note = modifier_matches(prefix_tokens, geo.title)
            scores.append(
                FieldScore(
                    "protein_modifier",
                    0.0,
                    s_score,
                    "pass" if s_score >= 1.0 else "fail",
                    observed=obs_sec,
                    expected="from_geo_title",
                    note=s_note,
                )
            )
        else:
            s_score, s_note = protein_matches(obs_sec, exp_sec)
            scores.append(
                FieldScore(
                    "protein_secondary",
                    0.0,
                    s_score,
                    "pass" if s_score >= 1.0 else "fail",
                    observed=obs_sec,
                    expected=exp_sec,
                    note=s_note,
                )
            )
            # Adjust primary protein score weight contribution for reCLIP dual-target names
            if s_score < 1.0 and p_score >= 1.0:
                scores[0].note += ";secondary_issue"

    # Cell line — authoritative: GEO characteristics (cell line:) > Flow source.
    # GSM title cell tokens are often wrong in this cohort; do not up-score against them.
    auth_cell = geo_parsed["geo_cell_authoritative"] or row.get("source", "")
    c_score, c_note = cell_matches(parsed["cell"], auth_cell)
    title_cell = geo_parsed["geo_cell_title"]
    if c_score < 1.0 and title_cell and cell_matches(parsed["cell"], title_cell)[0] >= 1.0:
        c_note = f"matches_geo_title_only:{title_cell};authoritative={auth_cell}"
    scores.append(
        FieldScore(
            "cell",
            FIELD_WEIGHTS["cell"],
            c_score,
            "pass" if c_score >= 1.0 else ("warn" if c_score >= 0.5 else "fail"),
            observed=parsed["cell"],
            expected=auth_cell,
            note=c_note,
        )
    )

    # Method
    flow_method = row.get("experimental_method", "")
    exp_method = geo_parsed["geo_method"] or flow_method
    m_score, m_note = method_matches(parsed["method"], exp_method)
    scores.append(
        FieldScore(
            "method",
            FIELD_WEIGHTS["method"],
            m_score,
            "pass" if m_score >= 1.0 else "fail",
            observed=parsed["method"],
            expected=exp_method,
            note=m_note,
        )
    )

    # Condition / UV / EGF
    cond_score, cond_note = condition_matches(
        parsed["condition_token"],
        row.get("condition", ""),
        geo.title,
        geo.treatment_protocol,
    )
    scores.append(
        FieldScore(
            "condition",
            FIELD_WEIGHTS["condition"],
            cond_score,
            "pass" if cond_score >= 1.0 else ("warn" if cond_score >= 0.5 else "fail"),
            observed=parsed["condition_token"],
            expected=row.get("condition", "") or "from_geo_title",
            note=cond_note,
        )
    )

    # Replicate
    exp_rep = geo_parsed["geo_rep"]
    obs_rep = parsed["rep"]
    if exp_rep and obs_rep:
        r_score = 1.0 if exp_rep == obs_rep else 0.0
        r_note = "exact" if r_score else f"geo_rep={exp_rep}"
    elif not exp_rep and not obs_rep:
        r_score, r_note = 1.0, "n/a"
    else:
        r_score, r_note = 0.5, "partial"
    scores.append(
        FieldScore(
            "replicate",
            FIELD_WEIGHTS["replicate"],
            r_score,
            "pass" if r_score >= 1.0 else "fail",
            observed=obs_rep,
            expected=exp_rep,
            note=r_note,
        )
    )

    # Section (molecular-weight fraction; GEO uses "Section N")
    exp_sec = geo_parsed["geo_section"]
    obs_sec = parsed["section"]
    if exp_sec and obs_sec:
        s_score = 1.0 if exp_sec == obs_sec else 0.0
        s_note = "exact" if s_score else f"geo_section={exp_sec}"
    elif not exp_sec and not obs_sec:
        s_score, s_note = 1.0, "n/a"
    else:
        s_score, s_note = 0.5, "partial"
    scores.append(
        FieldScore(
            "section",
            FIELD_WEIGHTS["section"],
            s_score,
            "pass" if s_score >= 1.0 else ("warn" if s_score >= 0.5 else "fail"),
            observed=obs_sec,
            expected=exp_sec,
            note=s_note,
        )
    )

    # HMW (high-molecular-weight) fraction flag in sample name
    obs_hmw = parsed.get("hmw", "")
    title_has_hmw = "hmw" in geo.title.lower()
    if obs_hmw or title_has_hmw:
        h_score = 1.0 if obs_hmw and title_has_hmw else 0.0
        scores.append(
            FieldScore(
                "hmw",
                0.0,
                h_score,
                "pass" if h_score >= 1.0 else "fail",
                observed=obs_hmw or "(none)",
                expected="HMW" if title_has_hmw else "(none)",
                note="exact" if h_score else "hmw_flag_mismatch",
            )
        )

    # Species
    sp_score = 1.0 if parsed["species"] == "Hs" and "homo sapiens" in geo.organism.lower() else 0.0
    scores.append(
        FieldScore(
            "species",
            FIELD_WEIGHTS["species"],
            sp_score,
            "pass" if sp_score >= 1.0 else "fail",
            observed=parsed["species"],
            expected=geo.organism,
            note="exact" if sp_score else "mismatch",
        )
    )

    weighted = [fs for fs in scores if fs.weight > 0]
    total_w = sum(fs.weight for fs in weighted)
    confidence = round(100.0 * sum(fs.score * fs.weight for fs in weighted) / total_w, 1) if total_w else 0.0

    fails = [fs for fs in scores if fs.weight > 0 and fs.score < 1.0]
    if any(fs.field == "cell" and fs.score == 0.0 for fs in fails):
        band = "LOW"
    elif confidence >= 95.0:
        band = "HIGH"
    elif confidence >= 80.0:
        band = "MEDIUM"
    else:
        band = "LOW"
    return scores, confidence, band


def audit_rows(rows: Sequence[Dict[str, str]], delay_s: float) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    cache: Dict[str, GeoGsm] = {}
    out_rows: List[Dict[str, Any]] = []
    for i, row in enumerate(rows):
        gsm = (row.get("geo") or "").strip()
        geo = fetch_gsm(gsm, cache, delay_s=delay_s if i else 0.0)
        parsed = parse_sample_name(row.get("sample_name", ""))
        geo_parsed = parse_geo_title(geo)
        field_scores, confidence, band = score_sample(row, geo, parsed, geo_parsed)

        issues = [
            f"{fs.field}:{fs.note}"
            for fs in field_scores
            if fs.weight > 0 and fs.score < 1.0
        ]
        secondary_issues = [
            f"{fs.field}:{fs.note}"
            for fs in field_scores
            if fs.weight == 0 and fs.score < 1.0
        ]

        out_rows.append(
            {
                "sample_id": row.get("sample_id", ""),
                "sample_name": row.get("sample_name", ""),
                "gsm": gsm,
                "geo_title": geo.title,
                "geo_cell_line": geo.cell_line,
                "geo_cell_title": geo_parsed["geo_cell_title"],
                "flow_source": row.get("source", ""),
                "flow_purification_target": row.get("purification_target", ""),
                "flow_condition": row.get("condition", ""),
                "flow_method": row.get("experimental_method", ""),
                "name_protein_primary": parsed.get("protein_primary", ""),
                "name_protein_secondary": parsed.get("protein_secondary", ""),
                "name_cell": parsed.get("cell", ""),
                "name_method": parsed.get("method", ""),
                "name_condition": parsed.get("condition_token", ""),
                "name_hmw": parsed.get("hmw", ""),
                "name_rep": parsed.get("rep", ""),
                "name_section": parsed.get("section", ""),
                "confidence_score": confidence,
                "confidence_band": band,
                "weighted_issues": "; ".join(issues),
                "secondary_issues": "; ".join(secondary_issues),
                "geo_title_vs_cell_line_conflict": (
                    "yes"
                    if geo_parsed["geo_cell_title"]
                    and geo.cell_line
                    and geo_parsed["geo_cell_title"] != geo.cell_line
                    else "no"
                ),
                "field_scores_json": json.dumps(
                    [{k: getattr(fs, k) for k in ("field", "weight", "score", "status", "observed", "expected", "note")} for fs in field_scores]
                ),
            }
        )
        if (i + 1) % 25 == 0:
            logging.info("Audited %d / %d", i + 1, len(rows))

    summary = {
        "sample_count": len(out_rows),
        "confidence_bands": {},
        "avg_confidence": round(sum(r["confidence_score"] for r in out_rows) / len(out_rows), 1) if out_rows else 0,
        "geo_title_cell_conflicts": sum(1 for r in out_rows if r["geo_title_vs_cell_line_conflict"] == "yes"),
        "low_confidence_samples": [r["sample_name"] for r in out_rows if r["confidence_band"] == "LOW"],
        "issue_counts": {},
    }
    for r in out_rows:
        summary["confidence_bands"][r["confidence_band"]] = summary["confidence_bands"].get(r["confidence_band"], 0) + 1
        for issue in filter(None, r["weighted_issues"].split("; ")):
            key = issue.split(":", 1)[0]
            summary["issue_counts"][key] = summary["issue_counts"].get(key, 0) + 1
    return out_rows, summary


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input-csv", required=True, help="Flow pull CSV (flow_public_samples_pull_v3 output)")
    ap.add_argument("--output-csv", required=True, help="Audit output CSV")
    ap.add_argument("--summary-json", default="", help="Optional summary JSON path")
    ap.add_argument("--gsm-delay", type=float, default=0.12, help="Delay between GSM fetches (seconds)")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")

    with open(args.input_csv, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        logging.error("No rows in %s", args.input_csv)
        return 2

    out_rows, summary = audit_rows(rows, delay_s=args.gsm_delay)
    fieldnames = list(out_rows[0].keys())
    with open(args.output_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(out_rows)

    summary_path = args.summary_json or str(Path(args.output_csv).with_suffix(".summary.json"))
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    logging.info(
        "Wrote %s (%d rows, avg confidence %.1f)",
        args.output_csv,
        len(out_rows),
        summary["avg_confidence"],
    )
    logging.info("Summary: %s", summary_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

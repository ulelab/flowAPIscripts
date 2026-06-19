#!/usr/bin/env python3
"""
Audit Flow.bio CLIP sample names against external metadata.

ID resolution priority:
  1. GSM (geo field)
  2. ENCFF (from file_names)
  3. ENA / ERR (ena field or file_names)

Skips RNA-Seq, Ribo-Seq, and samples without any resolvable ID.

Example:
  python3 flow_sample_name_audit_all.py \\
    --input-csv projects/flow_public_samples_pull_audit_all.csv \\
    --output-csv projects/flow_sample_name_audit_all.csv
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
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import requests

from flow_sample_name_audit import (
    GeoGsm,
    cell_matches,
    fetch_gsm,
    method_matches,
    norm_token,
    protein_matches,
)

ENCODE_BASE = "https://www.encodeproject.org"
ENA_FILEREPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"

GSM_WEIGHTS = {"cell": 25.0, "protein": 25.0, "condition": 20.0, "method": 15.0, "name_overlap": 15.0}
ENCFF_WEIGHTS = {"target": 30.0, "biosample": 25.0, "method": 20.0, "name_overlap": 25.0}
ENA_WEIGHTS = {"method": 25.0, "name_overlap": 40.0, "species": 15.0, "target_hint": 20.0}

ENCFF_RE = re.compile(r"(ENCFF[A-Z0-9]+)", re.I)
ENA_RE = re.compile(r"(ER[RZ]\d+)", re.I)
GSM_RE = re.compile(r"(GSM\d+)", re.I)

NOISE_TOKENS = {
    "hs", "mm", "gg", "rep", "clip", "iclip", "eclip", "seclip", "reclip", "input", "ip",
    "cells", "cell", "min", "egf", "uvc", "nouv", "section", "hmw", "fastq", "gz", "merged",
    "cleaned", "fixed", "single", "paired", "control", "sminput",
}


@dataclass
class FieldScore:
    field: str
    weight: float
    score: float
    note: str = ""


def is_excluded_assay(row: Dict[str, str]) -> bool:
    st = (row.get("sample_type") or "").strip().lower()
    em = (row.get("experimental_method") or "").strip().lower()
    if st in {"rna-seq", "rna seq", "ribo-seq", "riboseq"}:
        return True
    if "rna-seq" in st or "ribo" in st:
        return True
    if "rna-seq" in em or "ribo" in em:
        return True
    if "quantseq" in em.replace(" ", "").replace("-", ""):
        return True
    if row.get("ribosome_type") or row.get("ribosome_stabilisation_method"):
        return True
    return False


def is_clip_row(row: Dict[str, str]) -> bool:
    return (row.get("sample_type") or "").strip().upper() == "CLIP"


def extract_gsm(row: Dict[str, str]) -> str:
    geo = (row.get("geo") or "").strip()
    if geo.upper().startswith("GSM"):
        return geo.upper()
    m = GSM_RE.search(geo)
    return m.group(1).upper() if m else ""


def extract_encff(row: Dict[str, str]) -> str:
    m = ENCFF_RE.search(row.get("file_names") or "")
    return m.group(1).upper() if m else ""


def extract_ena(row: Dict[str, str]) -> str:
    for src in (row.get("ena") or "", row.get("file_names") or ""):
        m = ENA_RE.search(src)
        if m:
            return m.group(1).upper()
    return ""


def resolve_lookup(row: Dict[str, str]) -> Tuple[str, str]:
    gsm = extract_gsm(row)
    if gsm:
        return "gsm", gsm
    encff = extract_encff(row)
    if encff:
        return "encff", encff
    ena = extract_ena(row)
    if ena:
        return "ena", ena
    return "", ""


def token_overlap_score(name: str, *refs: str) -> Tuple[float, str]:
    name_tokens = {norm_token(t) for t in re.findall(r"[A-Za-z0-9]+", name) if len(t) > 1}
    name_tokens -= NOISE_TOKENS
    ref_text = " ".join(r for r in refs if r)
    ref_tokens = {norm_token(t) for t in re.findall(r"[A-Za-z0-9]+", ref_text) if len(t) > 1}
    ref_tokens -= NOISE_TOKENS
    if not name_tokens:
        return 0.5, "no_name_tokens"
    if not ref_tokens:
        return 0.5, "no_ref_tokens"
    overlap = len(name_tokens & ref_tokens) / len(name_tokens)
    return round(min(1.0, overlap * 1.15), 3), f"overlap={overlap:.2f}"


def parse_geo_catalog(title: str) -> str:
    m = re.match(r"^[^(]+?\s*\((\d+)\)\s+eCLIP", title or "", re.I)
    return m.group(1) if m else ""


def split_protein_catalog(protein_token: str, catalog_hint: str = "") -> Tuple[str, str]:
    if catalog_hint and protein_token.endswith(catalog_hint):
        return protein_token[: -len(catalog_hint)], catalog_hint
    m = re.match(r"^(?P<gene>[A-Za-z][A-Za-z0-9]*?)(?P<cat>\d{4,5})$", protein_token)
    if m:
        return m.group("gene"), m.group("cat")
    return protein_token, ""


def extract_protein_from_geo_title(title: str) -> str:
    if not title:
        return ""
    patterns = [
        r"^([^(,]+?)\s*\(\d+\)\s+eCLIP",
        r"^([^(,]+?)\s+iCLIP",
        r"^([^(,]+?)\s+irCLIPv2",
        r"^([^(,]+?)\s+Re-CLIP",
        r"eCLIP from ([A-Za-z0-9]+)",
        r"Control eCLIP from ([A-Za-z0-9]+)",
        r"^([A-Za-z0-9]+),",
    ]
    for pat in patterns:
        m = re.search(pat, title, re.I)
        if m:
            return m.group(1).strip()
    return title.split(",")[0].strip()


def method_matches_extended(observed: str, expected: str) -> Tuple[float, str]:
    o = norm_token(observed)
    e = norm_token(expected)
    if not o or not e:
        return 0.0, "missing"
    base_score, base_note = method_matches(observed, expected)
    if base_score >= 1.0:
        return base_score, base_note
    if o == e:
        return 1.0, "exact"
    if {"seclip", "eclip"} >= {o, e}:
        return 1.0, "seclip_eclip_equiv"
    if "irclip" in o and "irclip" in e:
        return 1.0, "irclip_family"
    if "iclip" in o and "iclip" in e:
        return 1.0, "iclip_family"
    return 0.0, "mismatch"


def extract_cell_from_geo_title(title: str) -> str:
    for pat in [r"in\s+(\w+)\s+cells", r"from\s+(\w+)\s*\(", r"from\s+(\w+)\s+cells"]:
        m = re.search(pat, title, re.I)
        if m:
            return m.group(1)
    return ""


def extract_condition_from_geo(title: str, treatment: str) -> str:
    t = f"{title} {treatment}".lower()
    if re.search(r"\binput\b", t):
        return "INPUT"
    if re.search(r"\bip\b", t):
        return "IP"
    if "nouv" in t or "no uv" in t:
        return "noUV"
    if "uv cross" in t or "uvc" in t or re.search(r"\buvc\b", t):
        return "UVC"
    for label, cond in [("0min egf", "0min EGF"), ("15min egf", "15min EGF"), ("30min egf", "30min EGF"), ("60min egf", "60min EGF")]:
        if label in t:
            return cond
    return ""


def extract_method_from_geo(title: str, description: str) -> str:
    blob = f"{title} {description}".lower()
    for method in ("re-clip", "irclipv2", "irclip", "seclip", "eclip", "iclip2", "iclip", "miclip", "par-clip", "clip"):
        if method.replace("-", "") in blob.replace("-", ""):
            if method == "re-clip":
                return "Re-CLIP"
            if method == "irclipv2":
                return "irCLIP"
            if method == "iclip2":
                return "iCLIP2"
            return method.upper() if method in {"iclip", "eclip", "seclip", "miclip"} else method
    return ""


def library_name_from_geo(geo: GeoGsm) -> str:
    for line in (geo.description or "").split("\n"):
        if "library name:" in line.lower():
            return line.split(":", 1)[-1].strip()
    if geo.description and "library name:" in geo.description.lower():
        m = re.search(r"library name:\s*(.+)", geo.description, re.I)
        if m:
            return m.group(1).strip()
    return ""


def weighted_confidence(scores: List[FieldScore]) -> float:
    w = sum(s.weight for s in scores)
    if not w:
        return 0.0
    return round(100.0 * sum(s.score * s.weight for s in scores) / w, 1)


def confidence_band(score: float, hard_fail: bool = False) -> str:
    if hard_fail:
        return "LOW"
    if score >= 95:
        return "HIGH"
    if score >= 80:
        return "MEDIUM"
    return "LOW"


def score_gsm(row: Dict[str, str], geo: GeoGsm) -> Tuple[List[FieldScore], List[str]]:
    issues: List[str] = []
    title = geo.title
    lib = library_name_from_geo(geo)
    catalog_hint = parse_geo_catalog(title)

    auth_cell = geo.cell_line or geo.source_name or extract_cell_from_geo_title(title)
    flow_cell = row.get("source", "")
    c_score, c_note = cell_matches(flow_cell or _cell_from_name(row["sample_name"]), auth_cell)
    if c_score < 1.0 and auth_cell:
        nc = _cell_from_name(row["sample_name"])
        if nc:
            c2, _ = cell_matches(nc, auth_cell)
            c_score = max(c_score, c2)
    scores = [FieldScore("cell", GSM_WEIGHTS["cell"], c_score, c_note)]

    geo_prot = extract_protein_from_geo_title(title)
    name_tok = row["sample_name"].split("_")[0] if row.get("sample_name") else ""
    name_prot, _ = split_protein_catalog(name_tok, catalog_hint)
    flow_cond = row.get("condition", "") or _condition_from_name(row["sample_name"])
    flow_prot = row.get("purification_target", "")
    if norm_token(flow_cond) == "input" or norm_token(flow_prot) == "sminput":
        flow_prot = name_prot
    elif not flow_prot or norm_token(flow_prot) == "sminput":
        flow_prot = name_prot
    else:
        flow_prot = flow_prot or name_prot
    p_score, p_note = protein_matches(name_prot, geo_prot)
    if p_score < 1.0:
        p2, p2_note = protein_matches(flow_prot, geo_prot)
        if p2 > p_score:
            p_score, p_note = p2, p2_note
    if p_score < 1.0 and geo_prot and name_prot:
        ng, eg = norm_token(name_prot), norm_token(geo_prot)
        if eg.startswith(ng) or ng.startswith(eg):
            p_score, p_note = 1.0, "gene_prefix"
    scores.append(FieldScore("protein", GSM_WEIGHTS["protein"], p_score, p_note))

    geo_cond = extract_condition_from_geo(title, geo.treatment_protocol)
    if geo_cond and flow_cond:
        cond_score = 1.0 if norm_token(geo_cond) == norm_token(flow_cond) else 0.0
        scores.append(FieldScore("condition", GSM_WEIGHTS["condition"], cond_score, f"{flow_cond} vs {geo_cond}"))
    elif geo_cond or flow_cond:
        scores.append(FieldScore("condition", GSM_WEIGHTS["condition"], 0.5, "partial"))
    else:
        scores.append(FieldScore("condition", GSM_WEIGHTS["condition"], 0.5, "n/a"))

    geo_method = extract_method_from_geo(title, geo.description) or _method_from_name(row["sample_name"])
    flow_method = row.get("experimental_method", "") or _method_from_name(row["sample_name"])
    if geo_method and flow_method:
        m_score, m_note = method_matches_extended(flow_method, geo_method)
        scores.append(FieldScore("method", GSM_WEIGHTS["method"], m_score, m_note))
    else:
        scores.append(FieldScore("method", GSM_WEIGHTS["method"], 0.5, "n/a"))

    n_score, n_note = token_overlap_score(row["sample_name"], title, lib, geo_prot, auth_cell)
    scores.append(FieldScore("name_overlap", GSM_WEIGHTS["name_overlap"], n_score, n_note))

    if geo.cell_line and extract_cell_from_geo_title(title) and geo.cell_line != extract_cell_from_geo_title(title):
        issues.append("geo_title_cell_conflict")
    if geo.fetch_error:
        issues.append(f"gsm_fetch_error:{geo.fetch_error}")
    return scores, issues


def _cell_from_name(name: str) -> str:
    known = ["HEK293T", "HEK293", "K562", "HepG2", "A431", "HCT116", "mESC", "nESC"]
    for c in known:
        if re.search(rf"_{c}_", name, re.I) or name.upper().startswith(c.upper()):
            return c
    m = re.search(r"_Hs_([^_]+)_", name)
    return m.group(1) if m else ""


def _protein_from_name(name: str, catalog_hint: str = "") -> str:
    tok = name.split("_")[0] if name else ""
    return split_protein_catalog(tok, catalog_hint)[0]


def _method_from_name(name: str) -> str:
    for tok in name.split("_"):
        low = tok.lower()
        for method in ("irclipv2", "irclip", "seclip", "eclip", "iclip2", "iclip", "reclip", "miclip"):
            if method in low.replace("-", ""):
                if method == "irclipv2":
                    return "irCLIP"
                if method == "iclip2":
                    return "iCLIP2"
                if method == "reclip":
                    return "Re-CLIP"
                return method.upper() if method in {"iclip", "eclip", "seclip", "miclip"} else method
    return ""


def _condition_from_name(name: str) -> str:
    if re.search(r"_INPUT_", name, re.I):
        return "INPUT"
    if re.search(r"_IP_", name, re.I):
        return "IP"
    return ""


def fetch_encff(encff: str, cache: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    if encff in cache:
        return cache[encff]
    out: Dict[str, Any] = {"encff": encff, "error": ""}
    try:
        r = requests.get(f"{ENCODE_BASE}/files/{encff}/", headers={"Accept": "application/json"}, timeout=60)
        r.raise_for_status()
        f = r.json()
        out.update(
            {
                "paired_end": str(f.get("paired_end", "")),
                "assay_title": f.get("assay_title", ""),
                "output_type": f.get("output_type", ""),
                "target": (f.get("target") or {}).get("label", "") if isinstance(f.get("target"), dict) else "",
                "status": f.get("status", ""),
            }
        )
        ds = f.get("dataset", "")
        if isinstance(ds, str) and ds:
            er = requests.get(f"{ENCODE_BASE}{ds}?format=json", headers={"Accept": "application/json"}, timeout=60)
            if er.ok:
                exp = er.json()
                bio = exp.get("biosample_ontology") or {}
                out["biosample"] = bio.get("term_name", "") if isinstance(bio, dict) else ""
    except requests.RequestException as exc:
        out["error"] = str(exc)
    cache[encff] = out
    return out


def expected_clip_method_from_paired_end(pe: str) -> str:
    if pe == "1":
        return "seCLIP"
    if pe == "2":
        return "eCLIP"
    return ""


def score_encff(row: Dict[str, str], meta: Dict[str, Any]) -> Tuple[List[FieldScore], List[str]]:
    issues: List[str] = []
    if meta.get("error"):
        issues.append(f"encode_error:{meta['error']}")
    target = meta.get("target", "")
    flow_target = row.get("purification_target", "") or _protein_from_name(row["sample_name"])
    t_score, t_note = protein_matches(flow_target, target)
    scores = [FieldScore("target", ENCFF_WEIGHTS["target"], t_score, t_note)]

    biosample = meta.get("biosample", "")
    flow_source = row.get("source", "") or _cell_from_name(row["sample_name"])
    b_score, b_note = cell_matches(flow_source, biosample)
    scores.append(FieldScore("biosample", ENCFF_WEIGHTS["biosample"], b_score, b_note))

    exp_method = expected_clip_method_from_paired_end(str(meta.get("paired_end", "")))
    flow_method = row.get("experimental_method", "")
    if exp_method and flow_method:
        m_score, m_note = method_matches_extended(flow_method, exp_method)
        scores.append(FieldScore("method", ENCFF_WEIGHTS["method"], m_score, m_note))
    else:
        scores.append(FieldScore("method", ENCFF_WEIGHTS["method"], 0.5, "n/a"))

    n_score, n_note = token_overlap_score(row["sample_name"], target, biosample)
    scores.append(FieldScore("name_overlap", ENCFF_WEIGHTS["name_overlap"], n_score, n_note))
    if meta.get("output_type") not in ("reads", ""):
        issues.append(f"output_type={meta.get('output_type')}")
    return scores, issues


def fetch_ena_run(err: str, cache: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    if err in cache:
        return cache[err]
    out: Dict[str, Any] = {"err": err, "error": ""}
    try:
        r = requests.get(
            ENA_FILEREPORT,
            params={
                "accession": err,
                "result": "read_run",
                "format": "json",
                "fields": "run_accession,sample_accession,experiment_accession,library_layout,scientific_name",
            },
            timeout=60,
        )
        r.raise_for_status()
        rows = r.json()
        if rows:
            out.update(rows[0])
            sam = rows[0].get("sample_accession", "")
            if sam:
                r2 = requests.get(
                    ENA_FILEREPORT,
                    params={
                        "accession": sam,
                        "result": "sample",
                        "format": "json",
                        "fields": "sample_accession,sample_alias,scientific_name,description",
                    },
                    timeout=60,
                )
                if r2.ok and r2.json():
                    out["sample_alias"] = r2.json()[0].get("sample_alias", "")
                    out["sample_description"] = r2.json()[0].get("description", "")
    except requests.RequestException as exc:
        out["error"] = str(exc)
    cache[err] = out
    return out


def score_ena(row: Dict[str, str], meta: Dict[str, Any]) -> Tuple[List[FieldScore], List[str]]:
    issues: List[str] = []
    if meta.get("error"):
        issues.append(f"ena_error:{meta['error']}")
    layout = (meta.get("library_layout") or "").upper()
    exp_method = "seCLIP" if layout == "SINGLE" else ("eCLIP" if layout == "PAIRED" else "")
    flow_method = row.get("experimental_method", "")
    if exp_method and flow_method:
        m_score, m_note = method_matches_extended(flow_method, exp_method)
    else:
        m_score, m_note = 0.5, "n/a"
    scores = [FieldScore("method", ENA_WEIGHTS["method"], m_score, m_note)]

    alias = meta.get("sample_alias", "")
    n_score, n_note = token_overlap_score(row["sample_name"], alias, meta.get("sample_description", ""))
    scores.append(FieldScore("name_overlap", ENA_WEIGHTS["name_overlap"], n_score, n_note))

    org = (meta.get("scientific_name") or "").lower()
    sp_score = 1.0 if "homo sapiens" in org and "_Hs_" in row["sample_name"] else (0.5 if not org else 0.0)
    scores.append(FieldScore("species", ENA_WEIGHTS["species"], sp_score, org))

    target = row.get("purification_target", "") or _protein_from_name(row["sample_name"])
    desc = (meta.get("sample_description") or "").lower()
    if target and target.lower() in desc:
        th_score = 1.0
    elif target and norm_token(target) in norm_token(desc):
        th_score = 0.5
    else:
        th_score = 0.5 if not target else 0.0
    scores.append(FieldScore("target_hint", ENA_WEIGHTS["target_hint"], th_score, target))

    return scores, issues


def load_gsm_cache(path: Path) -> Dict[str, GeoGsm]:
    if not path.is_file():
        return {}
    data = json.loads(path.read_text(encoding="utf-8"))
    return {k: GeoGsm(**v) for k, v in data.items()}


def save_gsm_cache(path: Path, cache: Dict[str, GeoGsm]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({k: asdict(v) for k, v in cache.items()}, indent=0), encoding="utf-8")


def audit_all(
    rows: Sequence[Dict[str, str]],
    gsm_cache: Dict[str, GeoGsm],
    encff_cache: Dict[str, Dict[str, Any]],
    ena_cache: Dict[str, Dict[str, Any]],
    gsm_delay: float,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, str]], Dict[str, Any]]:
    audited: List[Dict[str, Any]] = []
    skipped: List[Dict[str, str]] = []
    excluded = 0
    not_clip = 0

    clip_rows = []
    for row in rows:
        if not is_clip_row(row):
            not_clip += 1
            continue
        if is_excluded_assay(row):
            excluded += 1
            continue
        clip_rows.append(row)

    for i, row in enumerate(clip_rows):
        mode, acc = resolve_lookup(row)
        if not mode:
            skipped.append({**row, "skip_reason": "no_gsm_encff_or_ena"})
            continue

        issues: List[str] = []
        ref_title = ""
        if mode == "gsm":
            geo = fetch_gsm(acc, gsm_cache, delay_s=gsm_delay if acc not in gsm_cache else 0.0)
            scores, issues = score_gsm(row, geo)
            ref_title = geo.title
        elif mode == "encff":
            meta = fetch_encff(acc, encff_cache)
            scores, issues = score_encff(row, meta)
            ref_title = f"{meta.get('target','')} {meta.get('biosample','')}"
        else:
            meta = fetch_ena_run(acc, ena_cache)
            scores, issues = score_ena(row, meta)
            ref_title = meta.get("sample_alias", "")

        conf = weighted_confidence(scores)
        hard_fail = any(
            s.weight > 0 and s.score == 0.0
            for s in scores
            if s.field in {"cell", "target", "protein", "biosample"}
        )
        band = confidence_band(conf, hard_fail=hard_fail and conf < 70)
        weighted_issues = "; ".join(f"{s.field}:{s.note}" for s in scores if s.score < 1.0)

        audited.append(
            {
                "sample_id": row.get("sample_id", ""),
                "sample_name": row.get("sample_name", ""),
                "project_id": row.get("project_id", ""),
                "project_name": row.get("project_name", ""),
                "lookup_mode": mode,
                "lookup_accession": acc,
                "geo": row.get("geo", ""),
                "ena": row.get("ena", ""),
                "file_names": row.get("file_names", ""),
                "flow_source": row.get("source", ""),
                "flow_purification_target": row.get("purification_target", ""),
                "flow_condition": row.get("condition", ""),
                "flow_method": row.get("experimental_method", ""),
                "ref_title_or_summary": ref_title,
                "confidence_score": conf,
                "confidence_band": band,
                "weighted_issues": weighted_issues,
                "triage_flags": "; ".join(issues),
                "field_scores_json": json.dumps([s.__dict__ for s in scores]),
            }
        )
        if (i + 1) % 100 == 0:
            logging.info("Audited %d / %d CLIP rows (%d skipped so far)", i + 1, len(clip_rows), len(skipped))

    summary = {
        "input_rows": len(rows),
        "clip_rows": len(clip_rows),
        "audited": len(audited),
        "skipped_no_id": len(skipped),
        "excluded_rna_ribo": excluded,
        "not_clip_sample_type": not_clip,
        "avg_confidence": round(sum(r["confidence_score"] for r in audited) / len(audited), 1) if audited else 0,
        "confidence_bands": {},
        "lookup_modes": {},
        "low_confidence_samples": [r["sample_name"] for r in audited if r["confidence_band"] == "LOW"][:50],
    }
    for r in audited:
        summary["confidence_bands"][r["confidence_band"]] = summary["confidence_bands"].get(r["confidence_band"], 0) + 1
        summary["lookup_modes"][r["lookup_mode"]] = summary["lookup_modes"].get(r["lookup_mode"], 0) + 1
    return audited, skipped, summary


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input-csv", required=True)
    ap.add_argument("--output-csv", required=True)
    ap.add_argument("--skipped-csv", default="")
    ap.add_argument("--summary-json", default="")
    ap.add_argument("--gsm-cache", default="/home/mikej10/advbfx/projects/flow_gsm_cache.json")
    ap.add_argument("--gsm-delay", type=float, default=0.04)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")

    with open(args.input_csv, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    gsm_cache = load_gsm_cache(Path(args.gsm_cache))
    encff_cache: Dict[str, Dict[str, Any]] = {}
    ena_cache: Dict[str, Dict[str, Any]] = {}

    audited, skipped, summary = audit_all(rows, gsm_cache, encff_cache, ena_cache, args.gsm_delay)
    save_gsm_cache(Path(args.gsm_cache), gsm_cache)

    out = Path(args.output_csv)
    if audited:
        with out.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(audited[0].keys()))
            w.writeheader()
            w.writerows(audited)

    skipped_path = Path(args.skipped_csv or str(out.with_name(out.stem + "_skipped.csv")))
    if skipped:
        with skipped_path.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(skipped[0].keys()))
            w.writeheader()
            w.writerows(skipped)

    summary_path = Path(args.summary_json or str(out.with_suffix(".summary.json")))
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    logging.info("Wrote %s (%d audited, %d skipped)", out, len(audited), len(skipped))
    logging.info("Summary: %s", summary_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

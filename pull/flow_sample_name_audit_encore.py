#!/usr/bin/env python3
"""
Audit ENCORE / ENCODE-style Flow.bio sample names against GEO GSM metadata,
the Yeo ENCORE manifest, and SRA layout (single-end -> seCLIP).

ENCORE samples are GEO/SRA-based (SRR filenames). ENCFF IDs are resolved via
ENCODE API when a matching released file exists (dbxrefs SRA:SRR...); most
ENCORE-native runs will not have ENCFF records.

Example:
  python3 flow_sample_name_audit_encore.py \\
    --input-csv projects/ENCORE/project_422456844677891835_samples.csv \\
    --sra-table projects/ENCORE/SraRunTable.csv \\
    --yeo-manifest projects/ENCODE/yeo/All\\ eCLIP\\ data-Table\\ 1.tsv \\
    --output-csv projects/ENCORE/sample_name_audit_full.csv
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
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import requests

ENCODE_BASE = "https://www.encodeproject.org"

FIELD_WEIGHTS = {
    "protein": 22.0,
    "catalog_id": 8.0,
    "cell": 22.0,
    "condition": 18.0,
    "replicate": 10.0,
    "species": 5.0,
    "method": 15.0,
}

PROTEIN_ALIASES: Dict[str, List[str]] = {
    "ALYREF": ["ALYREF", "THOC4"],
    "GARS": ["GARS", "GARS1"],
    "RECQ1": ["RECQ1", "RECQL"],
    "FMRP": ["FMRP", "FMR1"],
    "TDP43": ["TDP43", "TARDBP"],
    "UBCH7": ["UBCH7", "UBE2L3"],
    "CACTIN": ["CACTIN", "ACTB"],
    "SRSF9": ["SRSF9", "SRSF9"],
    "FUBP3": ["FUBP3", "FUBP3"],
}


def norm_token(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", (s or "").strip().lower())


def protein_matches(observed: str, expected: str) -> Tuple[float, str]:
    o, e = norm_token(observed), norm_token(expected)
    if not o or not e:
        return 0.0, "missing"
    if o == e:
        return 1.0, "exact"
    for _canon, aliases in PROTEIN_ALIASES.items():
        alias_norm = {norm_token(a) for a in aliases}
        if o in alias_norm and e in alias_norm:
            return 1.0, "alias"
    return 0.0, "mismatch"


@dataclass
class GeoGsm:
    gsm: str
    title: str = ""
    cell_line: str = ""
    treatment: str = ""
    library_name: str = ""
    descriptions: List[str] = None  # type: ignore[assignment]
    organism: str = ""
    fetch_error: str = ""

    def __post_init__(self) -> None:
        if self.descriptions is None:
            self.descriptions = []


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
            elif line.startswith("!Sample_organism_ch1 = "):
                out.organism = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_description = "):
                out.descriptions.append(line.split("=", 1)[1].strip())
            elif line.startswith("!Sample_characteristics_ch1 = "):
                val = line.split("=", 1)[1].strip()
                if ":" in val:
                    key, v = val.split(":", 1)
                    key = key.strip().lower()
                    v = v.strip()
                    if key == "cell line":
                        out.cell_line = v
                    elif key == "treatment":
                        out.treatment = v
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        out.fetch_error = str(exc)
    for d in out.descriptions:
        if d.lower().startswith("library name:"):
            out.library_name = d.split(":", 1)[1].strip()
            break
    cache[gsm] = out
    if delay_s > 0:
        time.sleep(delay_s)
    return out


def parse_geo_title(title: str) -> Dict[str, str]:
    # "YTHDC1 (5096) eCLIP in K562 cells, IP replicate 1"
    m_prot = re.match(r"^([^(]+?)\s*\((\d+)\)\s+eCLIP", title, re.I)
    protein, catalog = "", ""
    if m_prot:
        protein = m_prot.group(1).strip()
        catalog = m_prot.group(2).strip()
    cell = ""
    m_cell = re.search(r"in\s+(\w+)\s+cells", title, re.I)
    if m_cell:
        cell = m_cell.group(1)
    condition = "INPUT" if re.search(r"\bINPUT\b", title, re.I) else ""
    if re.search(r"\bIP\b", title, re.I):
        condition = "IP"
    rep = ""
    m_rep = re.search(r"replicate\s+(\d+)", title, re.I)
    if m_rep:
        rep = m_rep.group(1)
    return {
        "geo_protein": protein,
        "geo_catalog": catalog,
        "geo_cell": cell,
        "geo_condition": condition,
        "geo_rep": rep,
    }


def parse_encore_sample_name(name: str) -> Dict[str, str]:
    m = re.match(
        r"^(?P<protein_token>[A-Za-z0-9]+)_Hs_(?P<cell>[^_]+)_(?P<condition>INPUT|IP)_rep(?P<rep>\d+)$",
        name,
    )
    if not m:
        return {"parse_error": name}
    return {
        "protein_token": m.group("protein_token"),
        "species": "Hs",
        "cell": m.group("cell"),
        "condition": m.group("condition"),
        "rep": m.group("rep"),
    }


def split_protein_catalog(protein_token: str, catalog_hint: str) -> Tuple[str, str]:
    if catalog_hint and protein_token.endswith(catalog_hint):
        return protein_token[: -len(catalog_hint)], catalog_hint
    m = re.match(r"^(?P<gene>[A-Za-z][A-Za-z0-9]*?)(?P<cat>\d{4,5})$", protein_token)
    if m:
        return m.group("gene"), m.group("cat")
    return protein_token, ""


def expected_method_from_layout(layout: str) -> str:
    if norm_token(layout) == "single":
        return "seCLIP"
    if norm_token(layout) in {"paired", "pair"}:
        return "eCLIP"
    return ""


def fetch_encode_file_by_srr(srr: str, cache: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    srr = srr.upper()
    if srr in cache:
        return cache[srr]
    out: Dict[str, Any] = {"srr": srr, "encff": "", "paired_end": "", "assay_title": "", "target": ""}
    # ENCODE search by accession filter on individual files isn't indexed by SRR;
    # use @id lookup via search on dbxrefs when available.
    try:
        resp = requests.get(
            f"{ENCODE_BASE}/search/",
            params={"type": "File", "file_format": "fastq", "format": "json", "limit": "all", "dbxrefs": f"SRA:{srr}"},
            headers={"Accept": "application/json"},
            timeout=60,
        )
        if resp.status_code == 200:
            graph = resp.json().get("@graph", [])
            if graph:
                f = graph[0]
                out.update(
                    {
                        "encff": f.get("accession", ""),
                        "paired_end": str(f.get("paired_end", "")),
                        "assay_title": f.get("assay_title", ""),
                        "target": (f.get("target") or {}).get("label", "") if isinstance(f.get("target"), dict) else "",
                    }
                )
    except requests.RequestException:
        pass
    cache[srr] = out
    return out


def load_sra_table(path: Path) -> Dict[str, Dict[str, str]]:
    by_run: Dict[str, Dict[str, str]] = {}
    with path.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            run = (row.get("Run") or "").strip().upper()
            if run:
                by_run[run] = row
    return by_run


def load_yeo_manifest(path: Path) -> Dict[str, Dict[str, str]]:
    by_key: Dict[str, Dict[str, str]] = {}
    if not path.exists():
        return by_key
    with path.open(encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            exp = (row.get("Experiment") or "").strip()
            if exp:
                by_key[exp] = row
    return by_key


def yeo_lookup(yeo: Dict[str, Dict[str, str]], protein: str, cell: str, catalog: str) -> Optional[Dict[str, str]]:
    candidates = [
        f"{protein}_{cell}_{catalog}",
        f"{protein}_{cell}_{catalog}".replace("GARS", "GARS1"),
    ]
    for key in candidates:
        if key in yeo:
            return yeo[key]
    for exp, row in yeo.items():
        if norm_token(row.get("Cells", "")) == norm_token(cell) and catalog and catalog in exp:
            if norm_token(row.get("RBP", "")) == norm_token(protein) or norm_token(row.get("RBP_official", "")) == norm_token(protein):
                return row
    return None


@dataclass
class FieldScore:
    field: str
    weight: float
    score: float
    status: str
    observed: str = ""
    expected: str = ""
    note: str = ""


def score_row(
    row: Dict[str, str],
    geo: GeoGsm,
    geo_parsed: Dict[str, str],
    parsed: Dict[str, str],
    sra_row: Dict[str, str],
    encode_meta: Dict[str, Any],
    yeo_row: Optional[Dict[str, str]],
) -> Tuple[List[FieldScore], float, str, List[str]]:
    triage: List[str] = []
    if parsed.get("parse_error"):
        return [FieldScore("parse", 100, 0, "fail", note=parsed["parse_error"])], 0.0, "PARSE_ERROR", ["parse_error"]

    protein_token = parsed["protein_token"]
    catalog_hint = geo_parsed.get("geo_catalog", "")
    name_protein, name_catalog = split_protein_catalog(protein_token, catalog_hint)
    geo_protein = geo_parsed.get("geo_protein", "")
    flow_target = row.get("purification_target", "")
    expected_protein = geo_protein
    if norm_token(parsed["condition"]) == "input":
        # INPUT rows: name still carries protein token; Flow target is SMInput
        expected_protein = geo_protein

    p_score, p_note = protein_matches(name_protein, expected_protein)
    if p_score < 1.0 and flow_target and norm_token(parsed["condition"]) == "ip":
        p_score2, p_note2 = protein_matches(name_protein, flow_target)
        if p_score2 > p_score:
            p_score, p_note = p_score2, f"flow_target:{p_note2}"
    scores = [
        FieldScore("protein", FIELD_WEIGHTS["protein"], p_score, "pass" if p_score >= 1.0 else "fail", name_protein, expected_protein, p_note)
    ]

    cat_score = 1.0 if name_catalog and catalog_hint and name_catalog == catalog_hint else (0.0 if catalog_hint else 0.5)
    cat_note = "exact" if cat_score == 1.0 else f"name={name_catalog},geo={catalog_hint}"
    scores.append(FieldScore("catalog_id", FIELD_WEIGHTS["catalog_id"], cat_score, "pass" if cat_score >= 1.0 else "fail", name_catalog, catalog_hint, cat_note))

    auth_cell = geo.cell_line or row.get("source", "")
    c_score, c_note = (1.0, "exact") if norm_token(parsed["cell"]) == norm_token(auth_cell) else (0.0, "mismatch")
    scores.append(FieldScore("cell", FIELD_WEIGHTS["cell"], c_score, "pass" if c_score >= 1.0 else "fail", parsed["cell"], auth_cell, c_note))

    geo_cond = geo_parsed.get("geo_condition", "")
    flow_cond = (row.get("condition") or "").strip().upper()
    name_cond = parsed["condition"]
    exp_cond = geo_cond or flow_cond
    cond_score = 1.0 if norm_token(name_cond) == norm_token(exp_cond) else 0.0
    if not flow_cond and geo_cond:
        triage.append("missing_flow_condition")
    scores.append(FieldScore("condition", FIELD_WEIGHTS["condition"], cond_score, "pass" if cond_score >= 1.0 else "fail", name_cond, exp_cond, "exact" if cond_score else "mismatch"))

    geo_rep = geo_parsed.get("geo_rep", "")
    rep_score = 1.0 if geo_rep and parsed["rep"] == geo_rep else (0.5 if not geo_rep else 0.0)
    scores.append(FieldScore("replicate", FIELD_WEIGHTS["replicate"], rep_score, "pass" if rep_score >= 1.0 else "fail", parsed["rep"], geo_rep, "exact" if rep_score == 1.0 else f"geo_rep={geo_rep}"))

    sp_score = 1.0 if parsed["species"] == "Hs" and "homo sapiens" in geo.organism.lower() else 0.0
    scores.append(FieldScore("species", FIELD_WEIGHTS["species"], sp_score, "pass" if sp_score >= 1.0 else "fail", parsed["species"], geo.organism, "exact" if sp_score else "mismatch"))

    layout = sra_row.get("LibraryLayout", "")
    expected_method = expected_method_from_layout(layout)
    flow_method = row.get("experimental_method", "")
    m_score = 1.0 if norm_token(flow_method) == norm_token(expected_method) else 0.0
    if norm_token(expected_method) == "seclip" and norm_token(flow_method) == "eclip":
        triage.append("method_should_be_seCLIP")
    scores.append(
        FieldScore("method", FIELD_WEIGHTS["method"], m_score, "pass" if m_score >= 1.0 else "fail", flow_method, expected_method, f"layout={layout}")
    )

    if yeo_row:
        official = yeo_row.get("RBP_official", "")
        if official and norm_token(official) != norm_token(name_protein) and protein_matches(name_protein, official)[0] < 1.0:
            triage.append(f"yeo_rbp_official={official}")

    if geo.library_name:
        canonical = geo.library_name  # e.g. YTHDC1_K562_5096_IP_rep1
        flow_style = f"{name_protein}{name_catalog}_{parsed['cell']}_{parsed['condition']}_rep{parsed['rep']}"
        if norm_token(canonical.replace("_", "")) != norm_token(flow_style.replace("_", "")):
            triage.append(f"library_name_diff:{canonical}")

    if encode_meta.get("encff"):
        triage.append(f"encode_file={encode_meta['encff']}")
    else:
        triage.append("no_encode_encff_for_srr")

    total_w = sum(s.weight for s in scores)
    confidence = round(100.0 * sum(s.score * s.weight for s in scores) / total_w, 1)
    if any(s.field == "cell" and s.score == 0 for s in scores):
        band = "LOW"
    elif confidence >= 95:
        band = "HIGH"
    elif confidence >= 80:
        band = "MEDIUM"
    else:
        band = "LOW"
    return scores, confidence, band, triage


def srr_from_row(row: Dict[str, str]) -> str:
    fn = row.get("file_names", "")
    m = re.search(r"(SRR\d+)", fn, re.I)
    return m.group(1).upper() if m else ""


def audit_rows(
    rows: Sequence[Dict[str, str]],
    sra_by_run: Dict[str, Dict[str, str]],
    yeo: Dict[str, Dict[str, str]],
    delay_s: float,
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    gsm_cache: Dict[str, GeoGsm] = {}
    enc_cache: Dict[str, Dict[str, Any]] = {}
    out: List[Dict[str, Any]] = []

    for i, row in enumerate(rows):
        gsm = (row.get("geo") or "").strip()
        srr = srr_from_row(row)
        geo = fetch_gsm(gsm, gsm_cache, delay_s=delay_s if i else 0.0)
        geo_parsed = parse_geo_title(geo.title)
        parsed = parse_encore_sample_name(row.get("sample_name", ""))
        sra_row = sra_by_run.get(srr, {})
        encode_meta = fetch_encode_file_by_srr(srr, enc_cache) if srr else {}
        catalog_hint = geo_parsed.get("geo_catalog", "")
        name_protein, name_catalog = split_protein_catalog(parsed.get("protein_token", ""), catalog_hint)
        yeo_row = yeo_lookup(yeo, name_protein, parsed.get("cell", ""), name_catalog) if parsed.get("protein_token") else None

        scores, confidence, band, triage = score_row(row, geo, geo_parsed, parsed, sra_row, encode_meta, yeo_row)
        issues = "; ".join(f"{s.field}:{s.note}" for s in scores if s.score < 1.0)

        out.append(
            {
                "sample_id": row.get("sample_id", ""),
                "sample_name": row.get("sample_name", ""),
                "gsm": gsm,
                "srr": srr,
                "geo_title": geo.title,
                "geo_library_name": geo.library_name,
                "geo_cell_line": geo.cell_line,
                "flow_source": row.get("source", ""),
                "flow_purification_target": row.get("purification_target", ""),
                "flow_condition": row.get("condition", ""),
                "flow_method": row.get("experimental_method", ""),
                "sra_library_layout": sra_row.get("LibraryLayout", ""),
                "expected_method": expected_method_from_layout(sra_row.get("LibraryLayout", "")),
                "encode_encff": encode_meta.get("encff", ""),
                "encode_paired_end": encode_meta.get("paired_end", ""),
                "yeo_experiment": yeo_row.get("Experiment", "") if yeo_row else "",
                "yeo_rbp_official": yeo_row.get("RBP_official", "") if yeo_row else "",
                "name_protein": name_protein,
                "name_catalog_id": name_catalog,
                "name_cell": parsed.get("cell", ""),
                "name_condition": parsed.get("condition", ""),
                "name_rep": parsed.get("rep", ""),
                "confidence_score": confidence,
                "confidence_band": band,
                "weighted_issues": issues,
                "triage_flags": "; ".join(triage),
                "field_scores_json": json.dumps(
                    [{k: getattr(s, k) for k in ("field", "weight", "score", "status", "observed", "expected", "note")} for s in scores]
                ),
            }
        )
        if (i + 1) % 50 == 0:
            logging.info("Audited %d / %d", i + 1, len(rows))

    summary = {
        "sample_count": len(out),
        "avg_confidence": round(sum(r["confidence_score"] for r in out) / len(out), 1) if out else 0,
        "confidence_bands": {},
        "method_should_be_seclip": sum(1 for r in out if "method_should_be_seCLIP" in r["triage_flags"]),
        "missing_flow_condition": sum(1 for r in out if "missing_flow_condition" in r["triage_flags"]),
        "no_encode_encff": sum(1 for r in out if "no_encode_encff_for_srr" in r["triage_flags"]),
        "encode_encff_found": sum(1 for r in out if r["encode_encff"]),
        "low_confidence": [r["sample_name"] for r in out if r["confidence_band"] == "LOW"],
        "issue_field_counts": {},
    }
    for r in out:
        summary["confidence_bands"][r["confidence_band"]] = summary["confidence_bands"].get(r["confidence_band"], 0) + 1
        for issue in filter(None, r["weighted_issues"].split("; ")):
            key = issue.split(":", 1)[0]
            summary["issue_field_counts"][key] = summary["issue_field_counts"].get(key, 0) + 1
    return out, summary


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input-csv", required=True)
    ap.add_argument("--output-csv", required=True)
    ap.add_argument("--sra-table", default="/home/mikej10/advbfx/projects/ENCORE/SraRunTable.csv")
    ap.add_argument("--yeo-manifest", default="/home/mikej10/advbfx/projects/ENCODE/yeo/All eCLIP data-Table 1.tsv")
    ap.add_argument("--summary-json", default="")
    ap.add_argument("--gsm-delay", type=float, default=0.1)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")

    with open(args.input_csv, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    sra_by_run = load_sra_table(Path(args.sra_table))
    yeo = load_yeo_manifest(Path(args.yeo_manifest))

    out_rows, summary = audit_rows(rows, sra_by_run, yeo, delay_s=args.gsm_delay)
    fieldnames = list(out_rows[0].keys())
    with open(args.output_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(out_rows)

    summary_path = args.summary_json or str(Path(args.output_csv).with_suffix(".summary.json"))
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    logging.info("Wrote %s (%d rows, avg confidence %.1f)", args.output_csv, len(out_rows), summary["avg_confidence"])
    logging.info("seCLIP fixes needed: %d", summary["method_should_be_seclip"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

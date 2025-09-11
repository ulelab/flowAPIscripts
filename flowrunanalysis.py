#!/usr/bin/env python3
import argparse
import getpass
import logging
import sys
from typing import Dict, List, Tuple
import requests
import numpy as np
import fnmatch

# -------------------------
# CLI
# -------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Run Flow.bio CLIP analysis (fetch + client-side filter by sample name)")
    p.add_argument("--pid", "--PID", dest="project_id", required=True,
                   help="Flow.bio Project ID (string)")
    p.add_argument("--filter", nargs=2, metavar=("KEY", "VALUE"), default=None,
                   help='Metadata filter. Supported now: --filter sample_name "*A"')
    p.add_argument("-n", "--num-chunks", type=int, default=1,
                   help="Split selected samples into N executions using numpy.array_split (default: 1)")
    p.add_argument("--dry-run", action="store_true",
                   help="Resolve everything and print payloads without submitting")
    p.add_argument("--verbose", action="store_true",
                   help="Enable DEBUG logging")
    p.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"], default=None,
                   help="Explicit log level (overrides --verbose)")
    return p.parse_args()

# -------------------------
# Logging
# -------------------------
def setup_logging(args):
    level = logging.INFO
    if args.verbose:
        level = logging.DEBUG
    if args.log_level:
        level = getattr(logging, args.log_level)
    logging.basicConfig(level=level, format="%(asctime)s | %(levelname)-8s | %(message)s")

# -------------------------
# Fixed CLIP pipeline settings
# -------------------------
PIPELINE_CLIP = {
    "prep_execution_id": "602442133516131481",
    "pipeline_id":       "960154035051242353",
    "pipeline_version":  "1.6",
}

# -------------------------
# HTTP helpers
# -------------------------
API_BASE = "https://api.flow.bio"

def login(session: requests.Session) -> str:
    username = input("Enter your username: ")
    password = getpass.getpass("Enter your password: ")
    r = session.post(f"{API_BASE}/login", json={"username": username, "password": password}, timeout=30)
    _raise_for_status(r)
    data = r.json()
    if "token" not in data:
        raise RuntimeError("Invalid username or password (no token in response)")
    return data["token"]

def _raise_for_status(resp: requests.Response):
    try:
        resp.raise_for_status()
    except requests.HTTPError as e:
        msg = f"HTTP {resp.status_code} error: {resp.text}"
        raise requests.HTTPError(msg) from e

# -------------------------
# Flow.bio queries
# -------------------------
def fetch_all_project_samples(session: requests.Session, token: str, project_id: str, page_size: int = 100) -> List[Dict]:
    headers = {"Authorization": f"Bearer {token}"}
    page = 1
    collected: List[Dict] = []
    while True:
        logging.debug("Fetching project samples page=%d count=%d", page, page_size)
        r = session.get(
            f"{API_BASE}/projects/{project_id}/samples",
            params={"page": page, "count": page_size},
            headers=headers,
            timeout=30,
        )
        _raise_for_status(r)
        payload = r.json()
        samples = payload.get("samples", [])
        if not samples:
            break
        collected.extend(samples)
        if len(samples) < page_size:
            break
        page += 1
    logging.info("Fetched %d samples from project %s", len(collected), project_id)
    return collected

def resolve_pipeline_version_id(session: requests.Session, token: str, pipeline_id: str, version_label: str) -> str:
    headers = {"Authorization": f"Bearer {token}"}
    r = session.get(f"{API_BASE}/pipelines/{pipeline_id}", headers=headers, timeout=30)
    _raise_for_status(r)
    pipeline = r.json()
    match = next((v for v in pipeline.get("versions", []) if v.get("name") == version_label), None)
    if not match:
        raise RuntimeError(f"Version {version_label!r} not found on pipeline {pipeline_id}")
    return match["id"]

def fetch_prep_execution(session: requests.Session, token: str, prep_execution_id: str) -> Dict:
    headers = {"Authorization": f"Bearer {token}"}
    r = session.get(f"{API_BASE}/executions/{prep_execution_id}", headers=headers, timeout=30)
    _raise_for_status(r)
    return r.json()

# -------------------------
# File binding
# -------------------------
FILE_MAP = {
    "fasta": "Homo_sapiens.GRCh38.fasta",
    "gtf": "Homo_sapiens.GRCh38.109.gtf",
    "smrna_fasta": "Homo_sapiens.GRCh38.smrna.fasta",
    "fasta_fai": "Homo_sapiens.GRCh38.fasta.fai",
    "chrom_sizes": "Homo_sapiens.GRCh38.fasta.sizes",
    "target_genome_index": "star",
    "smrna_genome_index": "bowtie",
    "smrna_fasta_fai": "Homo_sapiens.GRCh38.smrna.fasta.fai",
    "smrna_chrom_sizes": "Homo_sapiens.GRCh38.smrna.fasta.sizes",
    "longest_transcript": "longest_transcript.txt",
    "longest_transcript_fai": "longest_transcript.fai",
    "longest_transcript_gtf": "longest_transcript.gtf",
    "filtered_gtf": "Homo_sapiens_filtered.gtf",
    "seg_gtf": "Homo_sapiens_seg.gtf",
    "seg_filt_gtf": "Homo_sapiens_filtered_seg.gtf",
    "seg_resolved_gtf": "Homo_sapiens_filtered_seg_genicOtherfalse.resolved.gtf",
    "seg_resolved_gtf_genic": "Homo_sapiens_filtered_seg_genicOthertrue.resolved.gtf",
    "regions_gtf": "Homo_sapiens_regions.gtf.gz",
    "regions_filt_gtf": "Homo_sapiens_filtered_regions.gtf.gz",
    "regions_resolved_gtf": "Homo_sapiens_filtered_regions_genicOtherfalse.resolved.gtf",
    "regions_resolved_gtf_genic": "Homo_sapiens_filtered_regions_genicOthertrue.resolved.gtf",
}

def build_data_params_from_execution(execution: Dict, file_map: Dict[str, str]) -> Dict[str, str]:
    inputs = list((execution.get("data_params") or {}).values())
    all_files = inputs[:]
    for proc_ex in execution.get("process_executions", []):
        all_files.extend(proc_ex.get("downstream_data", []))
    idx: Dict[str, Dict] = {}
    for f in all_files:
        fn = f.get("filename")
        if fn and fn not in idx:
            idx[fn] = f
    missing: List[Tuple[str, str]] = []
    data_params: Dict[str, str] = {}
    for param_name, expected_filename in file_map.items():
        match = idx.get(expected_filename)
        if not match:
            missing.append((param_name, expected_filename))
            continue
        data_params[param_name] = match["id"]
    if missing:
        msg = "\n".join([f"  - {p}: {fn}" for p, fn in missing])
        raise RuntimeError(f"Missing required reference files:\n{msg}")
    return data_params

# -------------------------
# Client-side filter (sample_name glob)
# -------------------------
def filter_by_sample_name(samples: List[Dict], glob_expr: str | None) -> List[Dict]:
    if not glob_expr:
        return samples
    matched = [s for s in samples if fnmatch.fnmatch((s.get("name") or ""), glob_expr)]
    logging.info("Filter sample_name=%r matched %d / %d samples", glob_expr, len(matched), len(samples))
    return matched

# -------------------------
# Main
# -------------------------
def main():
    args = parse_args()
    setup_logging(args)

    project_id = args.project_id

    # Fixed CLIP pipeline settings
    prep_execution_id = PIPELINE_CLIP["prep_execution_id"]
    pipeline_id       = PIPELINE_CLIP["pipeline_id"]
    pipeline_version  = PIPELINE_CLIP["pipeline_version"]

    # HTTP session & auth
    session = requests.Session()
    token = login(session)
    headers = {"Authorization": f"Bearer {token}"}

    # Resolve pipeline version ID
    version_id = resolve_pipeline_version_id(session, token, pipeline_id, pipeline_version)
    logging.debug("Resolved pipeline version id=%s for label=%s", version_id, pipeline_version)

    # Prep execution & reference files
    ex = fetch_prep_execution(session, token, prep_execution_id)
    fileset_id = (ex.get("fileset") or {}).get("id")
    if not fileset_id:
        raise RuntimeError("Prep execution has no fileset id")
    data_params = build_data_params_from_execution(ex, FILE_MAP)

    # Fetch all samples for the project (no TSV), then client-side filter by sample name glob
    project_samples = fetch_all_project_samples(session, token, project_id, page_size=100)

    # parse --filter
    name_glob = None
    if args.filter:
        key, value = args.filter
        if key.lower() != "sample_name":
            raise SystemExit("Only --filter sample_name \"<glob>\" is supported in this iteration.")
        name_glob = value

    selected = filter_by_sample_name(project_samples, name_glob)

    # Print filtered list exactly like your working snippet
    print(f"\nFiltered {len(selected)} samples:")
    for s in selected:
        print(f"  {s['id']}: {s['name']}")

    if not selected:
        raise SystemExit("No samples selected after applying filter.")

    # Split into N execution batches
    n_chunks = max(1, int(args.num_chunks))
    chunks = [list(chunk) for chunk in np.array_split(np.array(selected, dtype=object), n_chunks)]
    logging.info("Prepared %d execution batch(es)", len(chunks))

    # Build and (optionally) submit one execution per chunk
    run_urls = []
    for i, chunk in enumerate(chunks, start=1):
        rows = [{
            "sample": s["id"],
            "values": {
                "group": s.get("name", ""),
                "replicate": "1",
            }
        } for s in chunk]

        payload = {
            "params": {
                "move_umi_to_header": "true",
                "umi_header_format": "NNNNNNNNNN",
                "umi_separator": "_",
                "skip_umi_dedupe": "false",
                "crosslink_position": "start",
            },
            "data_params": data_params,
            "csv_params": {"samplesheet": {"rows": rows, "paired": "both"}},
            "retries": None,
            "nextflow_version": "24.04.2",
            "fileset": fileset_id,
            "resequence_samples": False,
        }

        if args.dry_run:
            logging.info("DRY RUN: Batch %d/%d: %d samples", i, len(chunks), len(chunk))
            for s in chunk[:10]:
                logging.info("  %s | id=%s", s.get("name", ""), s.get("id", ""))
            continue

        if i == 1:
            proceed = input("Submit? (y/n): ").strip().lower()
            if proceed != "y":
                logging.info("Aborted by user.")
                sys.exit(0)

        r = session.post(
            f"{API_BASE}/pipelines/versions/{version_id}/run",
            headers=headers,
            json=payload,
            timeout=60,
        )
        _raise_for_status(r)
        run_id = r.json().get("id")
        if not run_id:
            raise RuntimeError(f"Submission succeeded but no run id in response: {r.text}")
        url = f"https://app.flow.bio/executions/{run_id}"
        run_urls.append(url)
        print(url)

    if args.dry_run:
        logging.info("DRY RUN complete: %d batch(es) prepared.", len(chunks))

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.error("%s", e)
        sys.exit(1)

#!/usr/bin/env python3
import argparse
import logging
import sys
from typing import Dict, List, Tuple
import numpy as np
import re
from flowbio import Client
import requests
import getpass

# -------------------------
# CLI Usage - python3 ./flowrunanalysis_flowbio.py --pid ######### --filter sample_name 'STAU2_HepG2.*$' -n #ofbatches --start-batch 1 --end-batch 10
# -------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Run Flow.bio CLIP analysis using flowbio library (fetch + client-side filter by sample name)")
    p.add_argument("--pid", "--PID", dest="project_id", required=True,
                   help="Flow.bio Project ID (string)")
    p.add_argument("--filter", nargs=2, metavar=("KEY", "VALUE"), default=None,
                   help='Metadata filter. Supported now: --filter sample_name "<regex>"')
    p.add_argument("-n", "--num-chunks", type=int, default=1,
                   help="Split selected samples into N executions using numpy.array_split (default: 1)")
    p.add_argument("--start-batch", type=int, default=1,
                   help="Start execution from batch number (1-based, default: 1)")
    p.add_argument("--end-batch", type=int, default=None,
                   help="End execution at batch number (1-based, default: all batches)")
    return p.parse_args()

# -------------------------
# Logging
# -------------------------
def setup_logging():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")

# -------------------------
# CLIP pipeline settings
# -------------------------
PIPELINE_CLIP = {
    "prep_execution_id": "602442133516131481",
    "pipeline_name": "CLIP-Seq",
    "pipeline_version": "1.7",
    "nextflow_version": "24.04.2",
}

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

# -------------------------
# REST helpers for data discovery (prep execution and samples)
# -------------------------
API_BASE = "https://api.flow.bio"

def _raise_for_status(resp: requests.Response):
    try:
        resp.raise_for_status()
    except requests.HTTPError as e:
        msg = f"HTTP {resp.status_code} error: {resp.text}"
        raise requests.HTTPError(msg) from e

def rest_login(session: requests.Session) -> str:
    username = input("Enter your username: ")
    password = getpass.getpass("Enter your password: ")
    r = session.post(f"{API_BASE}/login", json={"username": username, "password": password}, timeout=30)
    _raise_for_status(r)
    data = r.json()
    if "token" not in data:
        raise RuntimeError("Invalid username or password (no token in response)")
    return data["token"]

def fetch_prep_execution(session: requests.Session, token: str, prep_execution_id: str) -> Dict:
    headers = {"Authorization": f"Bearer {token}"}
    r = session.get(f"{API_BASE}/executions/{prep_execution_id}", headers=headers, timeout=30)
    _raise_for_status(r)
    return r.json()

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

# -------------------------
# Client-side filter (sample_name regex)
# -------------------------
def filter_by_sample_name(samples: List[Dict], regex_expr: str | None) -> List[Dict]:
    if not regex_expr:
        return samples
    try:
        pattern = re.compile(regex_expr)
    except re.error as e:
        raise SystemExit(f"Invalid regex for sample_name: {e}")
    matched = [s for s in samples if pattern.search((s.get("name") or ""))]
    logging.info("Filter sample_name=%r matched %d / %d samples", regex_expr, len(matched), len(samples))
    return matched

# -------------------------
# Main
# -------------------------
def main():
    args = parse_args()
    setup_logging()

    project_id = args.project_id

    # CLIP pipeline settings
    prep_execution_id = PIPELINE_CLIP["prep_execution_id"]
    pipeline_name = PIPELINE_CLIP["pipeline_name"]
    pipeline_version = PIPELINE_CLIP["pipeline_version"]
    nextflow_version = PIPELINE_CLIP["nextflow_version"]

    # Initialize flowbio client (for submission only)
    client = Client()

    # REST session & auth for discovery endpoints
    session = requests.Session()
    token = rest_login(session)

    # Prep execution & reference files via REST
    ex = fetch_prep_execution(session, token, prep_execution_id)
    data_params = build_data_params_from_execution(ex, FILE_MAP)

    # Fetch all samples for the project via REST
    project_samples = fetch_all_project_samples(session, token, project_id, page_size=100)

    # Parse --filter
    name_regex = None
    if args.filter:
        key, value = args.filter
        if key.lower() != "sample_name":
            raise SystemExit("Only --filter sample_name \"<regex>\" is supported in this iteration.")
        name_regex = value

    selected = filter_by_sample_name(project_samples, name_regex)

    # Print filtered list
    print(f"\nFiltered {len(selected)} samples:")
    for s in selected:
        print(f"  {s['id']}: {s['name']}")

    if not selected:
        raise SystemExit("No samples selected after applying filter.")

    # Split into N execution batches
    n_chunks = max(1, int(args.num_chunks))
    chunks = [list(chunk) for chunk in np.array_split(np.array(selected, dtype=object), n_chunks)]
    logging.info("Prepared %d execution batch(es)", len(chunks))

    # Determine which batches to execute
    start_batch = max(1, args.start_batch)
    end_batch = args.end_batch if args.end_batch is not None else len(chunks)
    end_batch = min(end_batch, len(chunks))
    
    if start_batch > len(chunks):
        raise SystemExit(f"Start batch {start_batch} exceeds total batches {len(chunks)}")
    
    logging.info("Will execute batches %d to %d (out of %d total batches)", start_batch, end_batch, len(chunks))

    # Build and submit one execution per chunk in the specified range
    run_urls = []
    for i, chunk in enumerate(chunks, start=1):
        # Skip batches outside the specified range
        if i < start_batch or i > end_batch:
            logging.info("Skipping batch %d (not in range %d-%d)", i, start_batch, end_batch)
            continue

        # Build sample_params in the format expected by flowbio
        sample_params = {}
        for s in chunk:
            sample_params[s["id"]] = {
                "group": s.get("name", ""),
                "replicate": "1",
            }

        # Pipeline parameters
        params = {
            "move_umi_to_header": "false",
            "umi_header_format": "NNN",
            "umi_separator": "RX:Z:",
            "skip_umi_dedupe": "false",
            "crosslink_position": "start",
            "encode_eclip": "false",
        }

        if i == 1:
            proceed = input("Submit? (y/n): ").strip().lower()
            if proceed != "y":
                logging.info("Aborted by user.")
                sys.exit(0)

        try:
            # Use flowbio client to run pipeline
            execution = client.run_pipeline(
                name=pipeline_name,
                version=pipeline_version,
                nextflow_version=nextflow_version,
                params=params,
                data_params=data_params,
                sample_params={"samples": sample_params}
            )
            
            run_id = execution.get("id")
            if not run_id:
                raise RuntimeError(f"Submission succeeded but no run id in response: {execution}")
            
            url = f"https://app.flow.bio/executions/{run_id}"
            run_urls.append(url)
            print(f"Batch {i}: {url}")
            
        except Exception as e:
            logging.error("Failed to submit batch %d: %s", i, e)
            continue

    logging.info("Completed submission of %d batches", len(run_urls))

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.error("%s", e)
        sys.exit(1)

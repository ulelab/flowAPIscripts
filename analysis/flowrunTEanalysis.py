#!/usr/bin/env python3
import argparse
import logging
import sys
from typing import Dict, List, Tuple
import numpy as np
import re
import requests

FLOWBIO_USERNAME = "username"
FLOWBIO_PASSWORD = "password"

# -------------------------
# CLI Usage - python3 ./flowrunTEanalysis.py --pid ######### --filter sample_name 'STAU2_HepG2.*$' -n #ofbatches --start-batch 1 --end-batch 10
# -------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Run Flow.bio TE analysis using flowbio library (fetch + client-side filter by sample name)")
    p.add_argument("--pid", "--PID", dest="project_id", required=True,
                   help="Flow.bio Project ID (string)")
    p.add_argument("--filter", nargs=2, metavar=("KEY", "VALUE"), action="append", default=None,
                   help='Metadata filter. Supported: --filter sample_name "<regex>", --filter experimental_method "<text>"')
    p.add_argument("-n", "--num-chunks", type=int, default=1,
                   help="Split selected samples into N executions using numpy.array_split (default: 1)")
    p.add_argument("--start-batch", type=int, default=1,
                   help="Start execution from batch number (1-based, default: 1)")
    p.add_argument("--end-batch", type=int, default=None,
                   help="End execution at batch number (1-based, default: all batches)")
    p.add_argument("--limit", type=str, default=None,
                   help="Limit samples after filtering. Use 'N' for first N samples, or 'M-N' for range (e.g., '15-37'). Default: no limit")
    return p.parse_args()

# -------------------------
# Logging
# -------------------------
def setup_logging():
    pass

# -------------------------
# CLIP pipeline settings
# -------------------------
PIPELINE_CLIP = {
    "prep_execution_id": "583150835173081055",
    "pipeline_id": "860300013917014252",
    "pipeline_name": "hanalysis-clipseq",
    "pipeline_version": "test",
    "nextflow_version": "25.10.4",
}

# -------------------------
# Data parameters (file IDs)
# -------------------------
DATA_PARAMS = {
    "gtf": "865333367351887680",
    "fasta": "320468299270948664",
    "seg_gtf": "657490829677407079",
    "fasta_fai": "416915275728809491",
    "regions_gtf": "322316217358102943",
    "filtered_gtf": "166527636003294230",
    "genome_index": "626326653018779728",
    "regions_filt_gtf": "607452590642180296",
    "ncrna_chrom_sizes": "754069987244417321",
    "genome_chrom_sizes": "655637113164580324",
    "ncrna_genome_index": "660565400267158951",
    "regions_resolved_gtf": "484888182031316412",
    "representative_transcript": "166527636003294230",
    "representative_transcript_fai": "530453389518046265",
    "representative_transcript_gtf": "166527636003294230",
    "ncrna_fasta": "472564573697657643",
    "ncrna_fasta_fai": "896151603580010259",
    "telescope_gtf": "396036320917779302",
    "tetranscripts_gtf": "607147647259993520",
}

# -------------------------
# REST helpers for data discovery (prep execution and samples)
# -------------------------
API_BASE = "https://api.flow.bio"

def rest_login(session: requests.Session) -> str:
    r = session.post(f"{API_BASE}/login", json={"username": FLOWBIO_USERNAME, "password": FLOWBIO_PASSWORD}, timeout=30)
    r.raise_for_status()
    return r.json()["token"]

def resolve_pipeline_version_id(session: requests.Session, token: str, pipeline_id: str, version_label: str) -> str:
    headers = {"Authorization": f"Bearer {token}"}
    r = session.get(f"{API_BASE}/pipelines/{pipeline_id}", headers=headers, timeout=30)
    r.raise_for_status()
    pipeline = r.json()
    return next((v["id"] for v in pipeline.get("versions", []) if v.get("name") == version_label), None)

def fetch_prep_execution(session: requests.Session, token: str, prep_execution_id: str) -> Dict:
    headers = {"Authorization": f"Bearer {token}"}
    r = session.get(f"{API_BASE}/executions/{prep_execution_id}", headers=headers, timeout=30)
    r.raise_for_status()
    return r.json()

def fetch_all_project_samples(session: requests.Session, token: str, project_id: str, page_size: int = 100, fetch_details: bool = False) -> List[Dict]:
    headers = {"Authorization": f"Bearer {token}"}
    page = 1
    collected: List[Dict] = []
    while True:
        r = session.get(
            f"{API_BASE}/projects/{project_id}/samples",
            params={"page": page, "count": page_size},
            headers=headers,
            timeout=30,
        )
        r.raise_for_status()
        payload = r.json()
        samples = payload.get("samples", [])
        if not samples:
            break
        
        # If fetch_details is True, get full sample details including metadata
        if fetch_details:
            for sample in samples:
                sample_id = sample.get("id")
                if sample_id:
                    detail_r = session.get(
                        f"{API_BASE}/samples/{sample_id}",
                        headers=headers,
                        timeout=30,
                    )
                    detail_r.raise_for_status()
                    detail_data = detail_r.json()
                    # Merge detail data into sample
                    sample.update(detail_data)
        
        collected.extend(samples)
        if len(samples) < page_size:
            break
        page += 1
    return collected

# -------------------------
# Limit parsing
# -------------------------
def parse_limit(limit_str: str | None) -> Tuple[int | None, int | None]:
    """
    Parse limit string into start and end indices (0-based).
    Returns (start, end) where end is exclusive.
    Examples:
    - "14" -> (0, 14)  # first 14 samples
    - "15-37" -> (14, 37)  # samples 15-37 (1-based becomes 0-based)
    """
    if not limit_str:
        return None, None
    
    if '-' in limit_str:
        # Range format: M-N
        try:
            start_str, end_str = limit_str.split('-', 1)
            start = int(start_str.strip()) - 1  # Convert to 0-based
            end = int(end_str.strip())  # Keep as 1-based for end
            if start < 0 or end <= start:
                raise ValueError("Invalid range: start must be >= 1 and end must be > start")
            return start, end
        except ValueError as e:
            raise SystemExit(f"Invalid range format '{limit_str}': {e}. Use format like '15-37'")
    else:
        # Single number format: N (first N samples)
        try:
            count = int(limit_str.strip())
            if count <= 0:
                raise ValueError("Count must be > 0")
            return 0, count
        except ValueError as e:
            raise SystemExit(f"Invalid limit format '{limit_str}': {e}. Use format like '14' or '15-37'")

# -------------------------
# Client-side filters
# -------------------------
def filter_by_sample_name(samples: List[Dict], regex_expr: str | None) -> List[Dict]:
    if not regex_expr:
        return samples
    try:
        pattern = re.compile(regex_expr)
    except re.error as e:
        raise SystemExit(f"Invalid regex for sample_name: {e}")
    matched = [s for s in samples if pattern.search((s.get("name") or ""))]
    return matched

def filter_by_experimental_method(samples: List[Dict], search_text: str | None) -> List[Dict]:
    """Filter samples where experimental method field matches search_text (case-insensitive)"""
    if not search_text:
        return samples

    matched = []
    found_methods = set()  # For debugging
    
    for s in samples:
        exp_method = None
        metadata = s.get("metadata", {})
        if isinstance(metadata, dict):
            exp_obj = (
                metadata.get("experimentalMethod")
                or metadata.get("experimental_method")
                or metadata.get("Experimental Method")
            )
            if isinstance(exp_obj, dict):
                exp_method = exp_obj.get("value")
            elif isinstance(exp_obj, str):
                exp_method = exp_obj

        if not exp_method:
            exp_method = s.get("experimental_method") or s.get("experimentalMethod") or s.get("Experimental Method")
        
        if exp_method:
            found_methods.add(str(exp_method))
            if search_text.lower() == str(exp_method).lower():
                matched.append(s)
    
    # Debug output
    if not matched and found_methods:
        print(f"\nDEBUG: Found experimental methods in samples: {sorted(found_methods)}")
        print(f"DEBUG: Searching for: '{search_text}'")
    elif not found_methods:
        print(f"\nDEBUG: No experimental method found in any sample metadata")
        print(f"DEBUG: Sample metadata keys (first sample): {list(samples[0].keys()) if samples else 'no samples'}")
        if samples and isinstance(samples[0].get("metadata"), dict):
            print(f"DEBUG: Metadata keys (first sample): {list(samples[0].get('metadata', {}).keys())}")

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
    pipeline_id = PIPELINE_CLIP["pipeline_id"]
    pipeline_name = PIPELINE_CLIP["pipeline_name"]
    pipeline_version = PIPELINE_CLIP["pipeline_version"]
    nextflow_version = PIPELINE_CLIP["nextflow_version"]

    # REST session & auth
    session = requests.Session()
    token = rest_login(session)
    headers = {"Authorization": f"Bearer {token}"}

    # Resolve pipeline version ID
    version_id = resolve_pipeline_version_id(session, token, pipeline_id, pipeline_version)

    # Prep execution & reference files via REST
    ex = fetch_prep_execution(session, token, prep_execution_id)
    fileset_id = ex["fileset"]["id"]
    data_params = DATA_PARAMS.copy()

    # Parse --filter to determine if we need full sample details
    needs_details = False
    if args.filter:
        for filter_key, filter_value in args.filter:
            if filter_key.lower() == "experimental_method":
                needs_details = True  # This is in full metadata, need to fetch details
                break
    
    # Fetch all samples for the project via REST (with details if needed for filtering)
    project_samples = fetch_all_project_samples(session, token, project_id, page_size=100, fetch_details=needs_details)

    # Apply filters
    selected = project_samples
    if args.filter:
        for key, value in args.filter:
            if key.lower() == "sample_name":
                selected = filter_by_sample_name(selected, value)
            elif key.lower() == "experimental_method":
                selected = filter_by_experimental_method(selected, value)
            else:
                raise SystemExit(f"Unsupported filter key: {key}. Supported: sample_name, experimental_method")

    # Apply limit if specified
    if args.limit:
        start_idx, end_idx = parse_limit(args.limit)
        original_count = len(selected)
        
        if start_idx is not None and end_idx is not None:
            # Validate range
            if start_idx >= len(selected):
                raise SystemExit(f"Start index {start_idx + 1} exceeds available samples ({len(selected)})")
            if end_idx > len(selected):
                end_idx = len(selected)
            
            selected = selected[start_idx:end_idx]

    # Print filtered list
    print(f"\nFiltered {len(selected)} samples:")
    for s in selected:
        print(f"  {s['id']}: {s['name']}")

    if not selected:
        raise SystemExit("No samples selected after applying filter.")

    # Split into N execution batches
    n_chunks = max(1, int(args.num_chunks))
    chunks = [list(chunk) for chunk in np.array_split(np.array(selected, dtype=object), n_chunks)]

    # Determine which batches to execute
    start_batch = max(1, args.start_batch)
    end_batch = args.end_batch if args.end_batch is not None else len(chunks)
    end_batch = min(end_batch, len(chunks))
    
    if start_batch > len(chunks):
        raise SystemExit(f"Start batch {start_batch} exceeds total batches {len(chunks)}")

    # Build and submit one execution per chunk in the specified range
    run_urls = []
    for i, chunk in enumerate(chunks, start=1):
        # Skip batches outside the specified range
        if i < start_batch or i > end_batch:
            continue

        # Build rows in the format expected by the pipeline
        rows = [{
            "sample": s["id"],
            "values": {
            }
        } for s in chunk]
        
        # Pipeline parameters
        params = {
            "move_umi_to_header": "false",
            #"umi_header_format": "NNNNNNNNNN",
            "umi_separator": "rbc:",
            "skip_umi_dedupe": "false",
            "crosslink_position": "start",
            "encode_eclip": "true",
            "run_te": "true",
            "source": "fastq",
            #"star_params": "--outFilterMultimapNmax 100 --outFilterMultimapScoreRange 1 --outSAMattributes All --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterType BySJout --alignIntronMin 20 --alignIntronMax 1000000 --outFilterScoreMin 10 --alignEndsType Extend5pOfRead1 --twopassMode Basic --limitOutSJcollapsed 4000000",
        }

        # Build payload for REST API submission
        payload = {
            "params": params,
            "data_params": data_params,
            "csv_params": {"input": {"rows": rows}},
            "retries": None,
            "nextflow_version": nextflow_version,
            "fileset": fileset_id,
            "resequence_samples": False,
        }

        if i == 1:
            proceed = input("Submit? (y/n): ").strip().lower()
            if proceed != "y":
                sys.exit(0)

        r = session.post(
            f"{API_BASE}/pipelines/versions/{version_id}/run",
            headers=headers,
            json=payload,
            timeout=120,
        )
        r.raise_for_status()
        run_id = r.json()["id"]
        
        url = f"https://app.flow.bio/executions/{run_id}"
        run_urls.append(url)
        print(f"Batch {i}: {url}")

if __name__ == "__main__":
    main()

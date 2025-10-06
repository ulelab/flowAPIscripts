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
# CLI Usage - python3 ./flowRNAanalysis.py --pid ######### --filter sample_name "*A" -n #ofbatches
# -------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Run Flow.bio RNA-seq analysis (fetch + client-side filter by sample name)")
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
# CLIP pipeline settings
# -------------------------
PIPELINE_RNA = {
    "prep_execution_id": "538918850957626478",
    "pipeline_id":       "583494301973770088",
    "pipeline_version":  "3.12",
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
    # Core genome references
    "fasta": "Homo_sapiens.GRCh38.fasta",
    "gtf": "Homo_sapiens.GRCh38.109.gtf",

    # RNA-seq specific references and indices
    "gene_bed": "Homo_sapiens.GRCh38.109.bed",
    "transcript_fasta": "genome.transcripts.fa",
    "splicesites": "Homo_sapiens.GRCh38.109.splice_sites.txt",

    # Aligner / quant indices (directories)
    "star_index": "star",
    "rsem_index": "rsem",
    "hisat2_index": "hisat2",
    "salmon_index": "salmon",
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

    # RNA pipeline settings
    prep_execution_id = PIPELINE_RNA["prep_execution_id"]
    pipeline_id       = PIPELINE_RNA["pipeline_id"]
    pipeline_version  = PIPELINE_RNA["pipeline_version"]

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

    # Build and submit one execution per chunk
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
                # UMI
                "with_umi": "true",
                "umitools_extract_method": "regex",
                "umitools_bc_pattern": "^(?P<discard_1>.{4})(?P<umi_1>.{5})",
                "skip_umi_extract": "false",
                "umitools_dedup_stats": "false",
                "save_umi_intermeds": "false",

                # Annotation/grouping
                "gencode": "false",
                "gtf_extra_attributes": "gene_name",
                "gtf_group_features": "gene_id",
                "featurecounts_group_type": "gene_biotype",
                "featurecounts_feature_type": "exon",

                # Align/quant
                "aligner": "star_salmon",
                "pseudo_aligner": "",
                "bam_csi_index": "false",
                "star_ignore_sjdbgtf": "false",
                "stringtie_ignore_gtf": "false",
                "save_unaligned": "false",
                "save_align_intermeds": "false",
                "skip_markduplicates": "false",
                "skip_alignment": "false",
                "skip_pseudo_alignment": "false",

                # Trimming
                "trimmer": "trimgalore",
                "skip_trimming": "false",
                "save_trimmed": "false",

                # rRNA options
                "remove_ribo_rna": "false",
                "ribo_database_manifest": "./assets/rrna-db-defaults.txt",
                "save_non_ribo_reads": "false",

                # QC
                "deseq2_vst": "true",
                "skip_bigwig": "false",
                "skip_stringtie": "false",
                "skip_fastqc": "false",
                "skip_preseq": "true",
                "skip_qualimap": "false",
                "skip_rseqc": "false",
                "skip_biotype_qc": "false",
                "skip_deseq2_qc": "false",
                "skip_multiqc": "false",
                "skip_qc": "false",
            },
            # File/data bindings are provided via data_params (from FILE_MAP)
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

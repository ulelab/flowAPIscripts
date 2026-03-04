#!/usr/bin/env python3
"""
Flow.bio upload_sample API example.
Usage: python3 uploadsample_flowbio_v5.py -p PROJECT_ID -r ROWS [--base-dir DIR] file.xlsx
"""

import argparse
import os
import math
from typing import Dict, Any
from flowbio import Client
import pandas as pd

FLOWBIO_USERNAME = "username"
FLOWBIO_PASSWORD = "password"

def _get(row: Dict[str, Any], key: str) -> str:
    val = row.get(key)
    if val is None or (isinstance(val, float) and math.isnan(val)):
        return ""
    return str(val).strip()


def build_metadata(row: Dict[str, Any], project_id: int) -> Dict[str, Any]:
    """Build metadata dict for upload_sample from Excel row."""
    metadata = {
        "project": project_id,
        "sample_type": _get(row, "Type") or "CLIP",
        "organism": _get(row, "Organism"),
        "scientist": _get(row, "Scientist"),
        "pi": _get(row, "PI"),
        "organisation": _get(row, "Organisation"),
        "purification_agent": _get(row, "Purification Agent"),
        "experimental_method": _get(row, "Experimental Method"),
        "condition": _get(row, "Condition"),
        "sequencer": _get(row, "Sequencer"),
        "comments": _get(row, "Comments"),
        "five_prime_barcode_sequence": _get(row, "5' Barcode Sequence"),
        "three_prime_barcode_sequence": _get(row, "3' Barcode Sequence"),
        "three_prime_adapter_name": _get(row, "3' Adapter Name"),
        "three_prime_adapter_sequence": _get(row, "3' Adapter Sequence"),
        "read1_primer": _get(row, "Read 1 Primer"),
        "read2_primer": _get(row, "Read 2 Primer"),
        "rt_primer": _get(row, "RT Primer"),
        "umi_barcode_sequence": _get(row, "UMI Barcode Sequence"),
        "umi_separator": _get(row, "UMI Separator"),
        "geo": _get(row, "GEO ID"),
        "ena": _get(row, "ENA ID"),
        "pubmed": _get(row, "PubMed ID"),
        "source": _get(row, "Source") or _get(row, "Cell or Tissue"),
        "source_annotation": _get(row, "Source Text"),
        "purification_target": _get(row, "Protein (Purification Target)") or _get(row, "Purification Target"),
        "purification_target_annotation": _get(row, "Purification Target Annotation"),
        "strandedness": _get(row, "Strandedness (Required)") or _get(row, "Strandedness"),
        "rna_selection_method": _get(row, "RNA Selection Method"),
        "ribosome_type": _get(row, "Ribosome Type"),
        "size_selection": _get(row, "Size selection") or _get(row, "Size Selection"),
        "separation_method": _get(row, "Separation Method"),
        "ribosome_stabilisation_method": _get(row, "Ribosome stabilization method"),
    }
    return {k: v for k, v in metadata.items() if v}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="Excel file (.xlsx)")
    parser.add_argument("-p", "--project", type=int, required=True)
    parser.add_argument("-r", "--rows", required=True, help="Rows e.g. 1-3,5")
    parser.add_argument("--base-dir", default=".")
    args = parser.parse_args()

    df = pd.read_excel(args.file)
    df.columns = [str(c).strip() for c in df.columns]
    rows = df.to_dict(orient="records")

    # Parse rows (e.g. "1-3,5" -> [1,2,3,5])
    selected = set()
    for part in args.rows.split(","):
        part = part.strip()
        if "-" in part:
            a, b = map(int, part.split("-"))
            selected.update(range(max(1, a), min(len(rows), b) + 1))
        else:
            selected.add(int(part))
    selected = sorted(s for s in selected if 1 <= s <= len(rows))

    client = Client()
    client.login(FLOWBIO_USERNAME, FLOWBIO_PASSWORD)

    for idx in selected:
        row = rows[idx - 1]
        name = _get(row, "Sample Name")
        path1 = _get(row, "File 1") or _get(row, "File")
        path2 = _get(row, "File 2")

        path1 = os.path.abspath(os.path.join(args.base_dir, path1)) if not os.path.isabs(path1) else path1
        path2 = os.path.abspath(os.path.join(args.base_dir, path2)) if path2 and not os.path.isabs(path2) else path2

        metadata = build_metadata(row, args.project)
        print(f"Uploading {name}...")

        if path2:
            sample = client.upload_sample(name, path1, path2, chunk_size=1_000_000, retries=5, progress=True, metadata=metadata)
        else:
            sample = client.upload_sample(name, path1, chunk_size=1_000_000, retries=5, progress=True, metadata=metadata)

        sid = sample.get("id") if isinstance(sample, dict) else getattr(sample, "id", None)
        print(f"  -> ID: {sid}")


if __name__ == "__main__":
    main()

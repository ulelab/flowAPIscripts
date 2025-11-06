#!/usr/bin/env python3
"""
Simplified Flow.bio sample upload script using flowbio library.
Based on https://docs.api.flow.bio/flowbio/upload
"""

import argparse
import getpass
import json
import os
import sys
from typing import Dict, Any, Optional
import math

from flowbio import Client

# Work around pandas attempting to import numexpr (incompatible with NumPy 2.x on some envs)
os.environ.setdefault("PANDAS_NO_NUMEXPR", "1")
import pandas as pd

# Controlled vocabulary mappings (aligned with working script)
RNA_STRANDEDNESS_MAP: Dict[str, str] = {
    "reverse": "reverse",
    "forward": "forward",
    "unstranded": "unstranded",
    "auto": "auto",
    # common aliases
    "rf": "reverse",
    "fr": "forward",
}

RNA_SELECTION_MAP: Dict[str, str] = {
    "ribominus": "ribominus",
    "polya": "polya",
    "targeted": "targeted",
}

# Default settings
DEFAULT_PROJECT_ID = 229379887183914226
DEFAULT_CHUNK_SIZE = 1_000_000

def parse_rows(spec: str, nrows: int) -> list[int]:
    """Parse row specification like '1-10,15,22-24' into list of 1-based indices."""
    out = set()
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            try:
                a, b = int(a), int(b)
                a = max(1, a)
                b = min(nrows, b)
                if a <= b:
                    out.update(range(a, b + 1))
            except ValueError:
                pass
        else:
            try:
                i = int(part)
                if 1 <= i <= nrows:
                    out.add(i)
            except ValueError:
                pass
    return sorted(out)

def normalize_path(path_from_csv: str, base_dir: str) -> str:
    """Convert relative paths to absolute paths."""
    p = os.path.expanduser((path_from_csv or "").strip())
    if not os.path.isabs(p):
        p = os.path.join(base_dir, p)
    return os.path.abspath(p)

def _get_cell_str(row: Dict[str, Any], key: str) -> str:
    """Safely retrieve a cell value as trimmed string; treat NaN/None as empty."""
    val = row.get(key)
    if val is None:
        return ""
    if isinstance(val, float) and math.isnan(val):
        return ""
    return str(val).strip()

def build_metadata(row: Dict[str, Any]) -> Dict[str, Any]:
    """Build metadata dictionary from TSV row, filtering out empty values."""
    metadata = {}
    
    # Required fields
    org = _get_cell_str(row, "Organism")
    if org:
        metadata["organism"] = org
    
    # Optional standard fields (mapping TSV column names to flowbio field names)
    field_mappings = {
        "Sample Name": "name",
        "Cell or Tissue": "source", 
        "Source": "source",
        "PI": "pi",
        "Scientist": "scientist",
        "Organisation": "organisation",
        "Organization": "organisation",
        "Institute": "organisation",
        "Purification Agent": "purificationAgent",
        "Condition": "condition",
        "Sequencer": "sequencer",
        "Comments": "comments",
        "Notes": "comments",
        "GEO ID": "geo",
        "GEO": "geo",
        "ENA ID": "ena", 
        "ENA": "ena",
        "PubMed ID": "pubmed",
        "PMID": "pubmed",
    }
    
    # Add mapped fields
    for tsv_col, flowbio_field in field_mappings.items():
        value = _get_cell_str(row, tsv_col)
        if value:
            metadata[flowbio_field] = value
    
    # Barcode and adapter fields (explicit camelCase mapping per Flow API docs)
    barcode_field_map = {
        "5' Barcode Sequence": "fivePrimeBarcodeSequence",
        "3' Barcode Sequence": "threePrimeBarcodeSequence",
        "3' Adapter Name": "threePrimeAdapterName",
        "3' Adapter Sequence": "threePrimeAdapterSequence",
        "RT Primer": "rtPrimer",
        "Read 1 Primer": "read1Primer",
        "Read 2 Primer": "read2Primer",
        "UMI Barcode Sequence": "umiBarcodeSequence",
        "UMI Separator": "umiSeparator",
    }

    for src, dst in barcode_field_map.items():
        value = _get_cell_str(row, src)
        if value:
            metadata[dst] = value
    
    # RNA-seq specific fields
    strandedness = _get_cell_str(row, "Strandedness").lower()
    if strandedness:
        metadata["strandedness"] = RNA_STRANDEDNESS_MAP.get(strandedness, strandedness)
    rna_sel = _get_cell_str(row, "RNA Selection Method").lower()
    if rna_sel:
        metadata["rnaSelectionMethod"] = RNA_SELECTION_MAP.get(rna_sel, rna_sel)
    
    return metadata

def main():
    parser = argparse.ArgumentParser(
        description="Simplified Flow.bio sample upload using flowbio library (Excel or TSV input)"
    )
    parser.add_argument("file", help="Path to input file (.xlsx, .xls, .tsv, or .txt)")
    parser.add_argument("--sheet", default=0, help="Worksheet name or index for Excel files (default: 0)")
    parser.add_argument("--rows", required=True, help="1-based rows to process, e.g. '1-10,15,22-24'")
    parser.add_argument("--base-dir", default=".", help="Base directory for relative file paths")
    parser.add_argument("--project-id", type=int, default=DEFAULT_PROJECT_ID, 
                       help=f"Project ID (default: {DEFAULT_PROJECT_ID})")
    parser.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE,
                       help=f"Upload chunk size in bytes (default: {DEFAULT_CHUNK_SIZE:,})")
    parser.add_argument("--retries", type=int, default=5, help="Number of upload retries (default: 5)")
    
    args = parser.parse_args()

    # Detect file type and read accordingly
    file_ext = os.path.splitext(args.file)[1].lower()
    is_excel = file_ext in ['.xlsx', '.xls']
    
    try:
        if is_excel:
            # Read Excel file
            sheet = args.sheet
            if isinstance(sheet, str) and sheet.isdigit():
                sheet = int(sheet)
            df = pd.read_excel(args.file, sheet_name=sheet)
            file_type_desc = f"Excel sheet '{args.sheet}'"
        else:
            # Read TSV file
            df = pd.read_csv(args.file, sep='\t')
            file_type_desc = "TSV file"
        
        # Normalize column names by stripping whitespace
        df.columns = [str(c).strip() for c in df.columns]
        rows = df.to_dict(orient="records")
    except Exception as e:
        file_type = "Excel" if is_excel else "TSV"
        print(f"Failed to read {file_type} file: {e}", file=sys.stderr)
        sys.exit(1)

    # Parse row selection
    selected_rows = parse_rows(args.rows, len(rows))
    if not selected_rows:
        print("No valid rows selected.", file=sys.stderr)
        sys.exit(1)

    print(f"Processing {len(selected_rows)} rows from {file_type_desc}...")

    # Get authentication
    username = input("Enter your Flow.bio username: ").strip()
    password = getpass.getpass("Enter your Flow.bio password: ")

    # Initialize client and login
    client = Client()
    try:
        client.login(username, password)
    except Exception as e:
        print(f"Login failed: {e}", file=sys.stderr)
        sys.exit(1)

    print("Successfully logged in to Flow.bio")

    # Process each selected row
    successful_uploads = 0
    failed_uploads = 0

    for idx in selected_rows:
        row = rows[idx - 1]  # Convert to 0-based index
        
        sample_name = _get_cell_str(row, "Sample Name")
        file_path = _get_cell_str(row, "File")
        
        if not sample_name:
            print(f"Row {idx}: Missing sample name, skipping", file=sys.stderr)
            failed_uploads += 1
            continue
            
        if not file_path:
            print(f"Row {idx} ({sample_name}): Missing file path, skipping", file=sys.stderr)
            failed_uploads += 1
            continue

        # Normalize file path
        full_path = normalize_path(file_path, args.base_dir)
        if not os.path.exists(full_path):
            print(f"Row {idx} ({sample_name}): File not found: {full_path}", file=sys.stderr)
            failed_uploads += 1
            continue

        # Build metadata
        metadata = build_metadata(row)
        metadata["project"] = args.project_id

        print(f"\nRow {idx}: Uploading '{sample_name}' from {full_path}")
        print(f"Metadata: {json.dumps(metadata, indent=2)}")

        try:
            # Upload sample using flowbio library
            result = client.upload_sample(
                name=sample_name,
                path1=full_path,
                chunk_size=args.chunk_size,
                progress=True,
                retries=args.retries,
                metadata=metadata
            )
            
            # Extract sample ID from result
            sample_id = result.get("id") if isinstance(result, dict) else getattr(result, "id", None)
            print(f"✓ Successfully uploaded '{sample_name}' (ID: {sample_id})")
            successful_uploads += 1
            
        except Exception as e:
            print(f"✗ Failed to upload '{sample_name}': {e}", file=sys.stderr)
            failed_uploads += 1

    # Summary
    print(f"\n{'='*50}")
    print(f"Upload Summary:")
    print(f"  Successful: {successful_uploads}")
    print(f"  Failed: {failed_uploads}")
    print(f"  Total: {successful_uploads + failed_uploads}")

if __name__ == "__main__":
    main()

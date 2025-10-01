#!/usr/bin/env python3
"""
Minimal Flow sample uploader (upload_sample only).

- Reads a TAB-delimited TSV.
- Select rows (1-based) via --rows like: 1-10,15,22-24.
- For each selected row, calls upload_sample with:
    name (from TSV), path1 (file), metadata {...} including many optional fields.
- Shows a progress bar and retries chunk uploads a few times internally.
- Treats a known GraphQL response bug as non-fatal when the upload actually succeeds.

Required TSV columns (fixed headers expected):
  Sample Name, File, Organism, 5' Barcode Sequence, Type
Optional TSV columns (if present, blanks are skipped):
  PI, Scientist, Organisation, Purification Agent, Experimental Method, Condition, Sequencer, Comments,
  3' Barcode Sequence, 3' Adapter Name, 3' Adapter Sequence, RT Primer, Read 1 Primer, Read 2 Primer,
  UMI Barcode Sequence, UMI Separator, Protein (Purification Target), GEO ID, ENA ID, PubMed ID, Cell or Tissue
"""

import argparse
import csv
import getpass
import os
import sys
from typing import List, Dict, Any

import flowbio  # available in your environment

PROJECT_ID = 931525950443849782  # your project
INTERNAL_RETRIES = 5             # fixed small retry count for chunk upload robustness
DEFAULT_CHUNK_SIZE = 1_000_000   # bytes


def parse_rows(spec: str, nrows: int) -> List[int]:
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


def norm_path(path_from_csv: str, base_dir: str) -> str:
    p = os.path.expanduser((path_from_csv or "").strip())
    if not os.path.isabs(p):
        p = os.path.join(base_dir, p)
    return os.path.abspath(p)


def add_if(meta: Dict[str, Any], key: str, val: str):
    val = (val or "").strip()
    if val:
        meta[key] = val


def main():
    ap = argparse.ArgumentParser(description="Minimal Flow upload_sample uploader from TSV")
    ap.add_argument("tsv", help="Tab-delimited annotation TSV")
    ap.add_argument("--rows", required=True, help="1-based rows, e.g. 1-10,15,22-24")
    ap.add_argument("--base-dir", default=".", help="Prepends to File when relative (default: .)")
    ap.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE, help=f"Chunk size bytes (default {DEFAULT_CHUNK_SIZE:,})")
    args = ap.parse_args()

    # Load TSV
    try:
        with open(args.tsv, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
    except Exception as e:
        print(f"[ERROR] Failed to read TSV: {e}", file=sys.stderr)
        sys.exit(2)
    if not rows:
        print("[ERROR] No rows in TSV.", file=sys.stderr)
        sys.exit(2)

    sel = parse_rows(args.rows, len(rows))
    if not sel:
        print("[ERROR] No valid rows selected.", file=sys.stderr)
        sys.exit(2)

    username = input("Enter your Flow.bio username: ").strip()
    password = getpass.getpass("Enter your Flow.bio password: ")

    client = flowbio.Client()
    try:
        client.login(username, password)
    except Exception as e:
        print(f"[ERROR] Login failed: {e}", file=sys.stderr)
        sys.exit(2)

    for idx in sel:
        row = rows[idx - 1]

        # Required fields
        sample_name = (row.get("Sample Name") or "").strip()
        file_col    = (row.get("File") or "").strip()
        organism    = (row.get("Organism") or "").strip()
        fivep       = (row.get("5' Barcode Sequence") or "").strip()
        category    = (row.get("Type") or "").strip()  # REQUIRED: populates 'category'

        if not sample_name:
            print(f"[WARN] Row {idx}: missing 'Sample Name'; skipping.", file=sys.stderr); continue
        if not file_col:
            print(f"[WARN] Row {idx} ({sample_name}): missing 'File'; skipping.", file=sys.stderr); continue
        if not organism:
            print(f"[WARN] Row {idx} ({sample_name}): missing 'Organism'; skipping.", file=sys.stderr); continue
        if not category:
            print(f"[WARN] Row {idx} ({sample_name}): missing 'Type' (category); skipping.", file=sys.stderr); continue
        if not fivep:
            # Five-prime barcode is required for CLIP; warn if absent.
            print(f"[WARN] Row {idx} ({sample_name}): missing \"5' Barcode Sequence\"; skipping.", file=sys.stderr); continue

        path1 = norm_path(file_col, args.base_dir)
        if not os.path.exists(path1):
            print(f"[WARN] Row {idx} ({sample_name}): file not found: {path1}; skipping.", file=sys.stderr); continue

        # Start with required fields (Flow expects camelCase keys)
        metadata: Dict[str, Any] = {
            "project": PROJECT_ID,          # numeric ID
            "organism": organism,
            "category": category,           # now REQUIRED and handled with other core fields
            "fivePrimeBarcodeSequence": fivep,
        }

        # Optional fields from fixed TSV headers
        add_if(metadata, "pi", row.get("PI", ""))
        add_if(metadata, "scientist", row.get("Scientist", ""))
        add_if(metadata, "organisation", row.get("Organisation", ""))

        # Purification agent (fills Purification Method box in UI)
        add_if(metadata, "purificationAgent", row.get("Purification Agent", ""))

        # Experimental method (adds detail to category)
        add_if(metadata, "experimentalMethod", row.get("Experimental Method", ""))

        add_if(metadata, "condition", row.get("Condition", ""))
        add_if(metadata, "sequencer", row.get("Sequencer", ""))
        add_if(metadata, "comments", row.get("Comments", ""))

        add_if(metadata, "threePrimeBarcodeSequence", row.get("3' Barcode Sequence", ""))
        add_if(metadata, "threePrimeAdapterName", row.get("3' Adapter Name", ""))
        add_if(metadata, "threePrimeAdapterSequence", row.get("3' Adapter Sequence", ""))
        add_if(metadata, "rtPrimer", row.get("RT Primer", ""))
        add_if(metadata, "read1Primer", row.get("Read 1 Primer", ""))
        add_if(metadata, "read2Primer", row.get("Read 2 Primer", ""))
        add_if(metadata, "umiBarcodeSequence", row.get("UMI Barcode Sequence", ""))
        add_if(metadata, "umiSeparator", row.get("UMI Separator", ""))

        # Cell or Tissue → source
        add_if(metadata, "source", row.get("Cell or Tissue", ""))

        # Protein target
        add_if(metadata, "purificationTarget", row.get("Protein (Purification Target)", ""))

        # IDs
        add_if(metadata, "geo", row.get("GEO ID", ""))
        add_if(metadata, "ena", row.get("ENA ID", ""))
        add_if(metadata, "pubmed", row.get("PubMed ID", ""))

        # Clean out blanks
        metadata = {k: v for k, v in metadata.items() if v not in ("", None)}

        try:
            s = client.upload_sample(
                sample_name,
                path1,
                progress=True,
                retries=INTERNAL_RETRIES,
                chunk_size=args.chunk_size,
                metadata=metadata,
            )
            sid = getattr(s, "id", None) or getattr(s, "sample_id", None)
            print(f"[OK] Uploaded row {idx}: '{sample_name}' (id={sid})")
        except Exception as e:
            emsg = str(e)
            benign = "Cannot query field 'type' on type 'SampleType'."
            if benign in emsg:
                print(f"[OK] Uploaded row {idx}: '{sample_name}' "
                      f"(server response quirk: {benign})")
                continue
            print(f"[ERROR] Upload failed for row {idx} '{sample_name}': {emsg}", file=sys.stderr)

    print("[DONE]")


if __name__ == "__main__":
    main()

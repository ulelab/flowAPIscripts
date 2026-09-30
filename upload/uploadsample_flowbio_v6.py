#!/usr/bin/env python3
"""
Upload sequencing samples to Flow using flowbio.v2 from an annotation sheet.

Supported formats: CSV, TSV, XLSX.

Example:
python3 uploadsample_flowbio_v6.py \
  --input test-datasets/annotation.csv \
  --rows 1-3 \
  --project-id 123 \
  --base-dir flowAPIscripts
"""

import argparse
import io
import math
import os
import sys
import warnings
from pathlib import Path
from typing import Any

os.environ.setdefault("PANDAS_NO_NUMEXPR", "1")

# Compatibility shim for Python versions without warnings.deprecated.
if not hasattr(warnings, "deprecated"):
    def _deprecated(_msg: str):  # type: ignore[no-redef]
        def _decorator(func):
            return func
        return _decorator
    warnings.deprecated = _deprecated  # type: ignore[attr-defined]

_saved_stderr = sys.stderr
try:
    sys.stderr = io.StringIO()
    import pandas as pd
finally:
    sys.stderr = _saved_stderr

from flowbio.v2 import Client, UsernamePasswordCredentials, ClientConfig


def _get(row: dict[str, Any], key: str) -> str:
    val = row.get(key)
    if val is None or (isinstance(val, float) and math.isnan(val)):
        return ""
    return str(val).strip()


def parse_rows(spec: str, nrows: int) -> list[int]:
    selected = set()
    for part in spec.split(","):
        token = part.strip()
        if not token:
            continue
        if "-" in token:
            start_s, end_s = token.split("-", 1)
            try:
                start_i, end_i = int(start_s), int(end_s)
            except ValueError:
                continue
            start_i = max(1, start_i)
            end_i = min(nrows, end_i)
            if start_i <= end_i:
                selected.update(range(start_i, end_i + 1))
        else:
            try:
                idx = int(token)
            except ValueError:
                continue
            if 1 <= idx <= nrows:
                selected.add(idx)
    return sorted(selected)


def resolve_path(path_str: str, base_dir: Path) -> Path:
    path_obj = Path(path_str).expanduser()
    if not path_obj.is_absolute():
        path_obj = base_dir / path_obj
    return path_obj.resolve()


def load_annotation_table(path: Path, sheet: str | int = 0) -> pd.DataFrame:
    """Load annotation rows from CSV, TSV, or XLSX."""
    suffix = path.suffix.lower()
    if suffix == ".xlsx":
        sheet_arg: str | int = int(sheet) if str(sheet).isdigit() else sheet
        df = pd.read_excel(path, sheet_name=sheet_arg)
    elif suffix == ".tsv":
        df = pd.read_csv(path, sep="\t")
    elif suffix in {".csv", ".txt"}:
        df = pd.read_csv(path)
    else:
        raise ValueError(
            f"Unsupported annotation format {path.suffix!r} for {path}. "
            "Use .csv, .tsv, or .xlsx"
        )
    df.columns = [str(c).strip() for c in df.columns]
    return df


def infer_data_files(row: dict[str, Any], base_dir: Path) -> dict[str, Path]:
    file1 = (
        _get(row, "File 1")
        or _get(row, "reads1")
        or _get(row, "Reads1")
        or _get(row, "File")
        or _get(row, "file")
    )
    file2 = _get(row, "File 2") or _get(row, "reads2") or _get(row, "Reads2")

    if not file1:
        raise ValueError("Missing reads1/file path (expected 'File 1' or 'File').")

    data = {"reads1": resolve_path(file1, base_dir)}
    if file2:
        data["reads2"] = resolve_path(file2, base_dir)
    return data


def collect_metadata(row: dict[str, Any]) -> dict[str, str]:
    metadata: dict[str, str] = {}
    for key, value in row.items():
        if value is None or (isinstance(value, float) and math.isnan(value)):
            continue
        key_s = str(key).strip()
        val_s = str(value).strip()
        if not key_s or not val_s:
            continue
        if key_s.startswith("metadata:"):
            metadata[key_s.split("metadata:", 1)[1].strip()] = val_s
        elif key_s.startswith("meta_"):
            metadata[key_s.split("meta_", 1)[1].strip()] = val_s

    # Pull metadata from known template columns (v5-compatible behavior).
    # Explicit metadata:* / meta_* keys above always win.
    column_to_metadata_key = {
        "Organism": "organism",
        "Scientist": "scientist",
        "PI": "pi",
        "Organisation": "organisation",
        "Purification Agent": "purification_agent",
        "Experimental Method": "experimental_method",
        "Condition": "condition",
        "Sequencer": "sequencer",
        "Comments": "comments",
        "5' Barcode Sequence": "five_prime_barcode_sequence",
        "3' Barcode Sequence": "three_prime_barcode_sequence",
        "3' Adapter Name": "three_prime_adapter_name",
        "3' Adapter Sequence": "three_prime_adapter_sequence",
        "Read 1 Primer": "read1_primer",
        "Read 2 Primer": "read2_primer",
        "RT Primer": "rt_primer",
        "UMI Barcode Sequence": "umi_barcode_sequence",
        "UMI Separator": "umi_separator",
        "GEO ID": "geo",
        "ENA ID": "ena",
        "PubMed ID": "pubmed",
        "Source": "source",
        "Cell or Tissue": "source",
        "Source Text": "source__annotation",
        "Protein (Purification Target)": "purification_target",
        "Purification Target": "purification_target",
        # Flow REST stores nested annotation on metadata objects as <key>__annotation
        # (see flow_public_samples_pull_v3.py / flow_public_samples_push_metadata_v2.py).
        "Purification Target Annotation": "purification_target__annotation",
        "purification_target__annotation": "purification_target__annotation",
        "source__annotation": "source__annotation",
        "Strandedness (Required)": "strandedness",
        "Strandedness": "strandedness",
        "RNA Selection Method": "rna_selection_method",
        "Ribosome Type": "ribosome_type",
        "Ribosome stabilization method": "ribosome_stabilisation_method",
        "Ribosome Stabilization Method": "ribosome_stabilisation_method",
        "Separation Method": "separation_method",
        "Size selection": "size_selection",
        "Size Selection": "size_selection",
    }
    for column_name, metadata_key in column_to_metadata_key.items():
        if metadata_key in metadata:
            continue
        value = _get(row, column_name)
        if value:
            metadata[metadata_key] = value

    return metadata


def main() -> int:
    parser = argparse.ArgumentParser(description="Flow.bio v2 sample upload from CSV/TSV/XLSX")
    parser.add_argument("--input", help="Path to annotation sheet (.csv, .tsv, or .xlsx)")
    parser.add_argument(
        "--input-xlsx",
        help="Deprecated alias for --input (XLSX/CSV/TSV all supported)",
    )
    parser.add_argument("--sheet", default=0, help="Sheet name or index for XLSX only (default: 0)")
    parser.add_argument("--rows", required=True, help="Row selection, e.g. 1-3,7")
    parser.add_argument("--project-id", required=True, help="Target Flow project ID")
    parser.add_argument("--base-dir", default=".", help="Base directory for relative file paths")
    parser.add_argument("--default-sample-type", default="RNA-Seq", help="Fallback sample type")
    parser.add_argument("--chunk-size", type=int, default=1_000_000, help="Upload chunk size bytes")
    parser.add_argument("--connection-retries", type=int, default=3, help="Connection retries")
    parser.add_argument("--show-progress", action="store_true", default=False, help="Show upload progress")
    parser.add_argument("--username", default=os.environ.get("FLOWBIO_USERNAME", ""), help="Flow username")
    parser.add_argument("--password", default=os.environ.get("FLOWBIO_PASSWORD", ""), help="Flow password")
    parser.add_argument("--dry-run", action="store_true", help="Validate inputs without uploading")
    args = parser.parse_args()

    input_path = args.input or args.input_xlsx
    if not input_path:
        print("Provide --input (or deprecated --input-xlsx).", file=sys.stderr)
        return 2

    if not args.username or not args.password:
        print("Missing credentials. Provide --username/--password or FLOWBIO_USERNAME/FLOWBIO_PASSWORD.", file=sys.stderr)
        return 2

    df = load_annotation_table(Path(input_path), sheet=args.sheet)
    rows = df.to_dict(orient="records")

    selected_rows = parse_rows(args.rows, len(rows))
    if not selected_rows:
        print("No valid rows selected.", file=sys.stderr)
        return 2

    base_dir = Path(args.base_dir).expanduser().resolve()
    config = ClientConfig(
        chunk_size=args.chunk_size,
        show_progress=args.show_progress,
        connection_retries=args.connection_retries,
    )
    client = Client(config=config)
    client.log_in(UsernamePasswordCredentials(username=args.username, password=args.password))

    success = 0
    failed = 0
    for idx in selected_rows:
        row = rows[idx - 1]
        sample_name = _get(row, "Sample Name") or _get(row, "Name") or f"row-{idx}"
        sample_type = _get(row, "Type") or _get(row, "sample_type") or args.default_sample_type

        try:
            data = infer_data_files(row, base_dir)
            for data_path in data.values():
                if not data_path.exists():
                    raise FileNotFoundError(f"Missing input file: {data_path}")

            metadata = collect_metadata(row)

            print(f"Row {idx}: {sample_name}")
            print(f"  sample_type={sample_type}")
            print(f"  data={{{', '.join([f'{k}: {v}' for k, v in data.items()])}}}")
            if metadata:
                print(f"  metadata keys={sorted(metadata.keys())}")

            if args.dry_run:
                success += 1
                continue

            sample = client.samples.upload_sample(
                name=sample_name,
                sample_type=sample_type,
                data=data,
                metadata=metadata or None,
                project_id=str(args.project_id),
            )
            print(f"  -> uploaded sample id={sample.id}")
            success += 1
        except Exception as exc:
            print(f"  -> failed: {exc}", file=sys.stderr)
            failed += 1

    print(f"\nCompleted. successful={success}, failed={failed}, total={len(selected_rows)}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""
Non-interactive Flow.bio sample upload script using flowbio library.
Credentials are embedded directly; use only in trusted environments.
"""

import argparse
import json
import os
import sys
from typing import Dict, Any, Optional, Set
import math

from flowbio import Client

# Prevent pandas from importing numexpr (can break with NumPy 2.x)
os.environ.setdefault("PANDAS_NO_NUMEXPR", "1")
import pandas as pd

# Controlled vocabulary mappings
RNA_STRANDEDNESS_MAP: Dict[str, str] = {
    "reverse": "reverse",
    "forward": "forward",
    "unstranded": "unstranded",
    "auto": "auto",
    "rf": "reverse",
    "fr": "forward",
}

RNA_SELECTION_MAP: Dict[str, str] = {
    "ribominus": "ribominus",
    "polya": "polya",
    "targeted": "targeted",
}

RIBO_TYPE_MAP: Dict[str, str] = {
    "monosome": "Monosome",
    "disome": "Disome",
    "polysomes": "Polysomes",
    "polysome": "Polysomes",
    "total ribosomes": "Total ribosomes",
    "totalribosomes": "Total ribosomes",
}

RIBO_SIZE_SELECTION_MAP: Dict[str, str] = {
    "standard monosome (28-30 nt)": "standard monosome (28-30 nt)",
    "broad selection (25-35 nt)": "broad selection (25-35 nt)",
    ">35 nt": "possible disomes (>35 nt)",
    "possible disomes (>35 nt)": "possible disomes (>35 nt)",
    "no size selection": "No size selection",
    "none": "No size selection",
    "other": "Other",
}

RIBO_SEPARATION_METHOD_MAP: Dict[str, str] = {
    "sucrose gradient": "Sucrose gradient",
    "cushion centrifugation": "Cushion centrifugation",
    "none": "None",
    "other": "Other",
}

UPDATE_SAMPLE_MUTATION = """
mutation UpdateSample(
  $id: ID!,
  $comments: String,
  $condition: String,
  $ena: String,
  $experimentalMethod: String,
  $fivePrimeBarcodeSequence: String,
  $geo: String,
  $name: String,
  $organisation: String,
  $organism: String,
  $pi: String,
  $private: Boolean,
  $project: String,
  $pubmed: String,
  $purificationAgent: String,
  $purificationTarget: String,
  $purificationTargetText: String,
  $read1Primer: String,
  $read2Primer: String,
  $rnaSelectionMethod: String,
  $rtPrimer: String,
  $scientist: String,
  $sequencer: String,
  $source: String,
  $sourceText: String,
  $strandedness: String,
  $threePrimeAdapterName: String,
  $threePrimeAdapterSequence: String,
  $threePrimeBarcodeSequence: String,
  $umiBarcodeSequence: String,
  $umiSeparator: String
) {
  updateSample(
    id: $id,
    comments: $comments,
    condition: $condition,
    ena: $ena,
    experimentalMethod: $experimentalMethod,
    fivePrimeBarcodeSequence: $fivePrimeBarcodeSequence,
    geo: $geo,
    name: $name,
    organisation: $organisation,
    organism: $organism,
    pi: $pi,
    private: $private,
    project: $project,
    pubmed: $pubmed,
    purificationAgent: $purificationAgent,
    purificationTarget: $purificationTarget,
    purificationTargetText: $purificationTargetText,
    read1Primer: $read1Primer,
    read2Primer: $read2Primer,
    rnaSelectionMethod: $rnaSelectionMethod,
    rtPrimer: $rtPrimer,
    scientist: $scientist,
    sequencer: $sequencer,
    source: $source,
    sourceText: $sourceText,
    strandedness: $strandedness,
    threePrimeAdapterName: $threePrimeAdapterName,
    threePrimeAdapterSequence: $threePrimeAdapterSequence,
    threePrimeBarcodeSequence: $threePrimeBarcodeSequence,
    umiBarcodeSequence: $umiBarcodeSequence,
    umiSeparator: $umiSeparator
  ) {
    sample { id }
  }
}
"""

# Default settings
DEFAULT_PROJECT_ID = 480752750995353723
DEFAULT_CHUNK_SIZE = 1_000_000
DEFAULT_SAMPLE_TYPE = "RNA-Seq"
ALLOWED_TSM: Dict[str, Set[str]] = {
    "CLIP": {
        "five_prime_barcode_sequence",
        "three_prime_barcode_sequence",
        "three_prime_adapter_name",
        "three_prime_adapter_sequence",
        "rt_primer",
        "read1_primer",
        "read2_primer",
        "umi_barcode_sequence",
        "umi_separator",
    },
    "RNA-SEQ": {
        "strandedness",
        "rna_selection_method",
        "three_prime_adapter_name",
        "three_prime_adapter_sequence",
        "read1_primer",
        "read2_primer",
    },
    "RIBO-SEQ": {
        "ribosome_type",
        "size_selection",
        "separation_method",
        "stabilization_method",
        "three_prime_adapter_name",
        "three_prime_adapter_sequence",
        "read1_primer",
        "read2_primer",
    },
}

# Static credentials (use with caution)
FLOWBIO_USERNAME = "username"
FLOWBIO_PASSWORD = "password"


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


def _normalize_vocab(value: str, mapping: Dict[str, str]) -> str:
    v = (value or "").strip()
    if not v:
        return ""
    if v.startswith("[") and v.endswith("]"):
        inner = v[1:-1].split(",", 1)[0]
        v = inner.strip().strip("'\"")
    v_lower = v.lower()
    if v_lower in mapping:
        return mapping[v_lower]

    compact = "".join(ch for ch in v_lower if ch.isalnum())
    for key, out in mapping.items():
        key_compact = "".join(ch for ch in key if ch.isalnum())
        if compact == key_compact:
            return out
    return v


def _build_type_specific_metadata(row: Dict[str, Any], sample_type: str) -> Dict[str, Any]:
    base: Dict[str, Any] = {
        "five_prime_barcode_sequence": _get_cell_str(row, "5' Barcode Sequence"),
        "three_prime_barcode_sequence": _get_cell_str(row, "3' Barcode Sequence"),
        "three_prime_adapter_name": _get_cell_str(row, "3' Adapter Name"),
        "three_prime_adapter_sequence": _get_cell_str(row, "3' Adapter Sequence"),
        "rt_primer": _get_cell_str(row, "RT Primer"),
        "read1_primer": _get_cell_str(row, "Read 1 Primer"),
        "read2_primer": _get_cell_str(row, "Read 2 Primer"),
        "umi_barcode_sequence": _get_cell_str(row, "UMI Barcode Sequence"),
        "umi_separator": _get_cell_str(row, "UMI Separator"),
        "strandedness": _normalize_vocab(_get_cell_str(row, "Strandedness"), RNA_STRANDEDNESS_MAP),
        "rna_selection_method": _normalize_vocab(_get_cell_str(row, "RNA Selection Method"), RNA_SELECTION_MAP),
        "ribosome_type": _normalize_vocab(_get_cell_str(row, "Ribosome Type"), RIBO_TYPE_MAP),
        "size_selection": _normalize_vocab(
            _get_cell_str(row, "Size selection") or _get_cell_str(row, "Size Selection"),
            RIBO_SIZE_SELECTION_MAP,
        ),
        "separation_method": _normalize_vocab(
            _get_cell_str(row, "Separation Method"),
            RIBO_SEPARATION_METHOD_MAP,
        ),
        "stabilization_method": _get_cell_str(row, "Ribosome stabilization method"),
    }

    tsm = {k: v for k, v in base.items() if v}
    st_key = (sample_type or "").strip().upper()
    allowed = ALLOWED_TSM.get(st_key)
    if allowed:
        tsm = {k: v for k, v in tsm.items() if k in allowed}
    else:
        safe_keys = {
            "three_prime_adapter_name",
            "three_prime_adapter_sequence",
            "read1_primer",
            "read2_primer",
        }
        tsm = {k: v for k, v in tsm.items() if k in safe_keys}
    return tsm


def _build_update_vars(sample_id: str, row: Dict[str, Any]) -> Dict[str, Any]:
    alias_map = {
        "name": ["Sample Name", "name"],
        "organism": ["Organism", "organism"],
        "source": ["Cell or Tissue", "Source", "source"],
        "experimentalMethod": ["Experimental Method", "experimental_method", "method"],
        "purificationAgent": ["Purification Agent", "purification_agent"],
        "purificationTarget": ["Protein (Purification Target)", "Purification Target", "purification_target"],
        "purificationTargetText": ["Purification Target Annotation", "purification_target_text"],
        "condition": ["Condition", "condition"],
        "sequencer": ["Sequencer", "sequencer", "platform"],
        "scientist": ["Scientist", "scientist"],
        "pi": ["PI", "pi"],
        "organisation": ["Organisation", "Organization", "Institute", "organisation"],
        "fivePrimeBarcodeSequence": ["5' Barcode Sequence", "five_prime_barcode_sequence"],
        "threePrimeBarcodeSequence": ["3' Barcode Sequence", "three_prime_barcode_sequence"],
        "threePrimeAdapterName": ["3' Adapter Name", "three_prime_adapter_name"],
        "threePrimeAdapterSequence": ["3' Adapter Sequence", "three_prime_adapter_sequence"],
        "read1Primer": ["Read 1 Primer", "read1_primer"],
        "read2Primer": ["Read 2 Primer", "read2_primer"],
        "rtPrimer": ["RT Primer", "rt_primer"],
        "umiBarcodeSequence": ["UMI Barcode Sequence", "umi_barcode_sequence"],
        "umiSeparator": ["UMI Separator", "umi_separator"],
        "strandedness": ["Strandedness", "strandedness"],
        "rnaSelectionMethod": ["RNA Selection Method", "rna_selection_method"],
        "geo": ["GEO ID", "GEO", "geo"],
        "ena": ["ENA ID", "ENA", "ena"],
        "pubmed": ["PubMed ID", "PMID", "pubmed"],
        "comments": ["Comments", "Notes", "comments"],
        "sourceText": ["Source Text"],
        "project": ["Project"],
        "private": ["private"],
    }

    def get_first(keys):
        for key in keys:
            val = _get_cell_str(row, key)
            if val:
                return val
        return ""

    vars_out: Dict[str, Any] = {"id": str(sample_id)}
    for target, keys in alias_map.items():
        val = get_first(keys)
        if not val:
            continue
        if target == "private":
            vars_out[target] = val.lower() in {"1", "true", "yes"}
        else:
            vars_out[target] = val

    vars_out.pop("project", None)
    return {k: v for k, v in vars_out.items() if v not in ("", None)}


def update_sample_metadata(client: Client, sample_id: str, row: Dict[str, Any]) -> None:
    if not hasattr(client, "execute"):
        return
    variables = _build_update_vars(sample_id, row)
    if len(variables) == 1:
        return
    try:
        client.execute(UPDATE_SAMPLE_MUTATION, variables=variables)
    except Exception as exc:
        print(f"[WARN] GraphQL updateSample failed for {sample_id}: {exc}", file=sys.stderr)


def build_metadata(row: Dict[str, Any], project_id: int) -> Dict[str, Any]:
    """Build upload metadata from TSV/Excel row for Flow.bio uploads."""
    sample_type = _get_cell_str(row, "Type") or DEFAULT_SAMPLE_TYPE

    metadata: Dict[str, Any] = {
        "project": project_id,
        "sample_type": sample_type,
    }

    org = _get_cell_str(row, "Organism")
    if org:
        metadata["organism"] = org

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

    for tsv_col, flowbio_field in field_mappings.items():
        value = _get_cell_str(row, tsv_col)
        if value:
            metadata[flowbio_field] = value

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

    exp_method = _get_cell_str(row, "Experimental Method")
    if exp_method:
        metadata["experimentalMethod"] = exp_method

    pur_target = _get_cell_str(row, "Protein (Purification Target)") or _get_cell_str(row, "Purification Target")
    if pur_target:
        metadata["purificationTarget"] = pur_target

    pur_target_text = _get_cell_str(row, "Purification Target Annotation")
    if pur_target_text:
        metadata["purificationTargetText"] = pur_target_text

    type_specific = _build_type_specific_metadata(row, sample_type)
    if type_specific:
        metadata["type_specific_metadata"] = json.dumps(type_specific)

    return {k: v for k, v in metadata.items() if v not in ("", None)}


def main():
    parser = argparse.ArgumentParser(
        description="Non-interactive Flow.bio sample upload (Excel or TSV input)"
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

    file_ext = os.path.splitext(args.file)[1].lower()
    is_excel = file_ext in ['.xlsx', '.xls']

    try:
        if is_excel:
            sheet = args.sheet
            if isinstance(sheet, str) and sheet.isdigit():
                sheet = int(sheet)
            df = pd.read_excel(args.file, sheet_name=sheet)
            file_type_desc = f"Excel sheet '{args.sheet}'"
        else:
            df = pd.read_csv(args.file, sep='\t')
            file_type_desc = "TSV file"

        df.columns = [str(c).strip() for c in df.columns]
        rows = df.to_dict(orient="records")
    except Exception as e:
        file_type = "Excel" if is_excel else "TSV"
        print(f"Failed to read {file_type} file: {e}", file=sys.stderr)
        sys.exit(1)

    selected_rows = parse_rows(args.rows, len(rows))
    if not selected_rows:
        print("No valid rows selected.", file=sys.stderr)
        sys.exit(1)

    print(f"Processing {len(selected_rows)} rows from {file_type_desc}...")

    client = Client()
    try:
        client.login(FLOWBIO_USERNAME, FLOWBIO_PASSWORD)
    except Exception as e:
        print(f"Login failed: {e}", file=sys.stderr)
        sys.exit(1)

    print("Successfully logged in to Flow.bio")
    if not hasattr(client, "execute"):
        for candidate in ("execute", "graphql", "gql"):
            func = getattr(client, candidate, None)
            if callable(func):
                setattr(client, "execute", func)
                if candidate != "execute":
                    print(f"[DEBUG] Using client.{candidate} for GraphQL execution.")
                break

    successful_uploads = 0
    failed_uploads = 0

    for idx in selected_rows:
        row = rows[idx - 1]

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

        full_path = normalize_path(file_path, args.base_dir)
        if not os.path.exists(full_path):
            print(f"Row {idx} ({sample_name}): File not found: {full_path}", file=sys.stderr)
            failed_uploads += 1
            continue

        metadata = build_metadata(row, args.project_id)

        print(f"\nRow {idx}: Uploading '{sample_name}' from {full_path}")
        print(f"Metadata: {json.dumps(metadata, indent=2)}")

        try:
            sample_flow = client.upload_sample(
                sample_name,
                full_path,
                chunk_size=args.chunk_size,
                retries=args.retries,
                progress=True,
                metadata=metadata,
            )

            sample_id = sample_flow.get("id") if isinstance(sample_flow, dict) else getattr(sample_flow, "id", None)
            if sample_id:
                update_sample_metadata(client, sample_id, row)
            print(f"✓ Successfully uploaded '{sample_name}' (ID: {sample_id})")
            successful_uploads += 1

        except Exception as e:
            print(f"✗ Failed to upload '{sample_name}': {e}", file=sys.stderr)
            failed_uploads += 1

    print(f"\n{'='*50}")
    print("Upload Summary:")
    print(f"  Successful: {successful_uploads}")
    print(f"  Failed: {failed_uploads}")
    print(f"  Total: {successful_uploads + failed_uploads}")


if __name__ == "__main__":
    main()


#!/usr/bin/env python3
"""
usage: uploadsample_flowbio_v4.py -p/--project PROJECT_ID -r/--rows ROWS [--base-dir BASE_DIR] annotationfile.xlsx
"""

import argparse
import json
import logging
import os
import sys
import io
import warnings
from contextlib import redirect_stdout
from typing import Dict, Any, Optional
import math
import requests
import re

from flowbio import Client

# Prevent pandas from importing numexpr (can break with NumPy 2.x)
os.environ.setdefault("PANDAS_NO_NUMEXPR", "1")
# Silence NumPy compatibility warnings
warnings.filterwarnings("ignore")
# Suppress stderr during pandas import to hide NumPy compatibility errors
_original_stderr = sys.stderr
_devnull = None
try:
    _devnull = open(os.devnull, 'w')
    sys.stderr = _devnull
    import pandas as pd
finally:
    if _devnull:
        _devnull.close()
    sys.stderr = _original_stderr

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

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


def normalize_path(path_from_excel: str, base_dir: str) -> str:
    """Convert relative paths to absolute paths."""
    p = os.path.expanduser((path_from_excel or "").strip())
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




def log_api_call(operation: str, request_data: Dict[str, Any], response: Any = None):
    """Log API calls for testing and debugging."""
    logger.info(f"\n{'='*60}")
    logger.info(f"API CALL: {operation}")
    logger.info(f"{'='*60}")
    logger.info(f"REQUEST:")
    logger.info(json.dumps(request_data, indent=2))
    if response is not None:
        logger.info(f"\nRESPONSE:")
        if isinstance(response, dict):
            logger.info(json.dumps(response, indent=2))
        else:
            logger.info(f"Response type: {type(response)}")
            logger.info(f"Response: {response}")
    logger.info(f"{'='*60}\n")


def build_upload_metadata(row: Dict[str, Any], project_id: int) -> Dict[str, Any]:
    # Get sample type from row
    sample_type = _get_cell_str(row, "Type")
    if not sample_type:
        sample_type = "CLIP"
        logger.warning("No 'Type' column found, defaulting to CLIP")
    
    # Start with required fields
    metadata: Dict[str, Any] = {
        "project": project_id,
        "sample_type": sample_type,
    }
    
    # Only send common fields in initial upload_sample() call
    basic_field_mappings = {
        "Sample Name": "name",
        "Organism": "organism",
        "PI": "pi",
        "Scientist": "scientist",
        "Organisation": "organisation",
        "Condition": "condition",
        "Sequencer": "sequencer",
        "Comments": "comments",
        "GEO ID": "geo",
        "ENA ID": "ena",
        "PubMed ID": "pubmed",
    }
    
    for tsv_col, flowbio_field in basic_field_mappings.items():
        value = _get_cell_str(row, tsv_col)
        if value:
            metadata[flowbio_field] = value
    
    source = _get_cell_str(row, "Source")
    if source:
        metadata["source"] = source
    
    source_text = _get_cell_str(row, "Source Text")
    if source_text:
        metadata["sourceText"] = source_text
    
    # Filter out empty values
    return {k: v for k, v in metadata.items() if v not in ("", None)}


def upload_sample_basic(client: Client, sample_name: str, file_path: str, 
                       metadata: Dict[str, Any]) -> Optional[str]:
    """
    Upload sample and extract sample_id from response.
    """
    logger.info(f"Uploading sample: {sample_name}")
    logger.info(f"File: {file_path}")
    
    log_api_call("upload_sample", {
        "sample_name": sample_name,
        "file_path": file_path,
        "metadata": metadata
    })
    
    # Capture stdout to extract sample_id from upload response
    captured_output = io.StringIO()
    
    try:
        with redirect_stdout(captured_output):
            client.upload_sample(
                sample_name,
                file_path,
                progress=True,
                metadata=metadata,
            )
    except Exception as e:
        # GraphQL error after successful upload - extract sample_id from stdout
        if "Cannot query field" not in str(e) and "GraphQlError" not in str(type(e)):
            logger.error(f"Upload failed: {e}", exc_info=True)
            log_api_call("upload_sample", metadata, {"error": str(e)})
            return None
    
    # Extract sample_id from stdout (works whether exception occurred or not)
    stdout_text = captured_output.getvalue()
    match = re.search(r'"sample_id"\s*:\s*"(\d+)"', stdout_text)
    if match:
        logger.info(f"✓ Extracted sample_id: {match.group(1)}")
        return match.group(1)
    
    logger.warning("Could not extract sample_id from stdout")
    return None


def update_sample_metadata_graphql(client: Client, sample_id: str, row: Dict[str, Any], 
                                   additional_metadata: Dict[str, Any] = None) -> bool:
    """
    Update sample metadata using GraphQL updateSample mutation.
    Returns True if successful, False otherwise.
    """
    if not hasattr(client, "execute"):
        logger.warning("Client does not have 'execute' method for GraphQL")
        return False
    
    # Combined field mapping: Excel columns -> Flow.bio field name -> GraphQL variable declaration
    # Format: {flowbio_field: {"excel_cols": [...], "graphql_var": "$fieldName: String"}}
    field_definitions = {
        "purificationAgent": {
            "excel_cols": ["Purification Agent"],
            "graphql_var": "$purificationAgent: String"
        },
        "experimentalMethod": {
            "excel_cols": ["Experimental Method"],
            "graphql_var": "$experimentalMethod: String"
        },
        "purificationTarget": {
            "excel_cols": ["Protein (Purification Target)"],
            "graphql_var": "$purificationTarget: String"
        },
        "purificationTargetText": {
            "excel_cols": ["Purification Target Annotation"],
            "graphql_var": "$purificationTargetText: String"
        },
        "fivePrimeBarcodeSequence": {
            "excel_cols": ["5' Barcode Sequence"],
            "graphql_var": "$fivePrimeBarcodeSequence: String"
        },
        "threePrimeBarcodeSequence": {
            "excel_cols": ["3' Barcode Sequence"],
            "graphql_var": "$threePrimeBarcodeSequence: String"
        },
        "threePrimeAdapterName": {
            "excel_cols": ["3' Adapter Name"],
            "graphql_var": "$threePrimeAdapterName: String"
        },
        "threePrimeAdapterSequence": {
            "excel_cols": ["3' Adapter Sequence"],
            "graphql_var": "$threePrimeAdapterSequence: String"
        },
        "rtPrimer": {
            "excel_cols": ["RT Primer"],
            "graphql_var": "$rtPrimer: String"
        },
        "read1Primer": {
            "excel_cols": ["Read 1 Primer"],
            "graphql_var": "$read1Primer: String"
        },
        "read2Primer": {
            "excel_cols": ["Read 2 Primer"],
            "graphql_var": "$read2Primer: String"
        },
        "umiBarcodeSequence": {
            "excel_cols": ["UMI Barcode Sequence"],
            "graphql_var": "$umiBarcodeSequence: String"
        },
        "umiSeparator": {
            "excel_cols": ["UMI Separator"],
            "graphql_var": "$umiSeparator: String"
        },
        "strandedness": {
            "excel_cols": ["Strandedness (Required)"],
            "graphql_var": "$strandedness: String"
        },
        "rnaSelectionMethod": {
            "excel_cols": ["RNA Selection Method"],
            "graphql_var": "$rnaSelectionMethod: String"
        },
        "ribosomeType": {
            "excel_cols": ["Ribosome Type"],
            "graphql_var": "$ribosomeType: String"
        },
        "sizeSelection": {
            "excel_cols": ["Size selection", "Size Selection"],
            "graphql_var": "$sizeSelection: String"
        },
        "separationMethod": {
            "excel_cols": ["Separation Method"],
            "graphql_var": "$separationMethod: String"
        },
        "ribosomeStabilisationMethod": {
            "excel_cols": ["Ribosome Stabilization Method"],
            "graphql_var": "$ribosomeStabilisationMethod: String"
        },
    }
    
    # Extract values from Excel row and build mutation variables
    mutation_variables = {"id": sample_id}
    fields_to_include = []
    mutation_args = []
    
    for flowbio_field, field_def in field_definitions.items():
        # Try to get value from Excel columns
        value = None
        for col in field_def["excel_cols"]:
            cell_value = _get_cell_str(row, col)
            if cell_value:
                value = cell_value
                break
        
        # If value found, add to mutation
        if value:
            mutation_variables[flowbio_field] = value
            fields_to_include.append(field_def["graphql_var"])
            mutation_args.append(f"{flowbio_field}: ${flowbio_field}")
    
    # Remove id from check (it's required, not optional)
    if len(mutation_variables) == 1:  # Only "id" present
        return True
    
    # Dynamically build GraphQL mutation with only the fields we have
    if not fields_to_include:
        return True
    
    var_declarations = ",\n      ".join(fields_to_include)
    mutation_args_str = ",\n        ".join(mutation_args)
    mutation = f"""
    mutation UpdateSample(
      $id: ID!,
      {var_declarations}
    ) {{
      updateSample(
        id: $id,
        {mutation_args_str}
      ) {{
        sample {{
          id
          name
        }}
      }}
    }}
    """
    
    try:
        result = client.execute(mutation, variables=mutation_variables)
        
        if result and result.get("data", {}).get("updateSample"):
            return True
        elif result and "errors" in result:
            return False
        else:
            return False
            
    except Exception:
        return False

def main():
    parser = argparse.ArgumentParser(
        description="Flow.bio sample upload script v4"
    )
    parser.add_argument("file", help="Path to input Excel file (.xlsx)")
    parser.add_argument("-p", "--project", type=int, required=True,
                       help="Project ID (required)")
    parser.add_argument("-r", "--rows", required=True,
                       help="1-based rows to process, e.g. '1-10,15,22-24'")
    parser.add_argument("--base-dir", default=".",
                       help="Base directory for relative file paths")
    parser.add_argument("--debug", action="store_true",
                       help="Enable debug logging")
    
    args = parser.parse_args()
    
    # Set logging level
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Validate file extension
    if not args.file.lower().endswith('.xlsx'):
        logger.error("Only .xlsx files are supported")
        sys.exit(1)
    
    # Read Excel file
    try:
        logger.info(f"Reading Excel file: {args.file}")
        df = pd.read_excel(args.file)
        df.columns = [str(c).strip() for c in df.columns]
        rows = df.to_dict(orient="records")
        
        logger.info(f"Loaded {len(rows)} rows")
        
    except Exception as e:
        logger.error(f"Failed to read Excel file: {e}", exc_info=True)
        sys.exit(1)
    
    # Parse row selection
    selected_rows = parse_rows(args.rows, len(rows))
    if not selected_rows:
        logger.error("No valid rows selected")
        sys.exit(1)
    
    logger.info(f"Processing {len(selected_rows)} row(s): {selected_rows}")
    
    # Initialize FlowBio client
    logger.info("Initializing FlowBio client...")
    client = Client()
    try:
        logger.info(f"Logging in with username: {FLOWBIO_USERNAME}")
        client.login(FLOWBIO_USERNAME, FLOWBIO_PASSWORD)
        logger.info("✓ Successfully logged in to Flow.bio")
    except Exception as e:
        logger.error(f"Login failed: {e}", exc_info=True)
        sys.exit(1)
    
    # Process selected rows
    successful = 0
    failed = 0
    
    for idx in selected_rows:
        row = rows[idx - 1]  # Convert to 0-based index
        
        logger.info(f"Processing row {idx}")
        
        # Extract basic info
        sample_name = _get_cell_str(row, "Sample Name")
        file_path = _get_cell_str(row, "File")
        
        logger.info(f"Sample Name: {sample_name}")
        logger.info(f"File: {file_path}")
        
        # Validate required fields
        if not sample_name:
            logger.error(f"Row {idx}: Missing sample name, skipping")
            failed += 1
            continue
        
        if not file_path:
            logger.error(f"Row {idx} ({sample_name}): Missing file path, skipping")
            failed += 1
            continue
        
        # Resolve file path
        full_path = normalize_path(file_path, args.base_dir)
        if not os.path.exists(full_path):
            logger.error(f"Row {idx} ({sample_name}): File not found: {full_path}")
            failed += 1
            continue
        
        # Get sample type from row
        sample_type = _get_cell_str(row, "Type")
        if sample_type:
            logger.info(f"Sample Type: {sample_type}")
        else:
            logger.warning(f"Row {idx} ({sample_name}): No 'Type' column found, defaulting to CLIP")
        
        # Build comprehensive metadata
        metadata = build_upload_metadata(row, args.project)
        
        # Upload sample
        sample_id = upload_sample_basic(
            client,
            sample_name,
            full_path,
            metadata
        )
        
        if sample_id:
            logger.info(f"✓ Successfully uploaded '{sample_name}'")
            logger.info(f"  Sample ID: {sample_id}")
            
            # Update missing metadata fields using GraphQL updateSample mutation
            logger.info(f"Updating metadata for sample {sample_id}...")
            update_success = update_sample_metadata_graphql(client, sample_id, row)
            
            if update_success:
                logger.info(f"✓ Successfully updated metadata for '{sample_name}'")
            else:
                logger.warning(f"⚠ Metadata update had issues for '{sample_name}'")
            
            successful += 1
        elif sample_id is None:
            # Upload succeeded but sample_id extraction failed
            # Try to extract from response text pattern
            logger.warning(f"Upload succeeded for '{sample_name}' but sample_id extraction failed")
            logger.warning(f"Check logs above for 'Response text' line to find sample_id")
            logger.warning(f"Metadata update skipped - sample_id required")
            successful += 1
        else:
            logger.error(f"✗ Failed to upload '{sample_name}'")
            failed += 1
    
    # Summary
    logger.info("Upload Summary:")
    logger.info(f"  Successful: {successful}")
    logger.info(f"  Failed: {failed}")
    logger.info(f"  Total: {successful + failed}")


if __name__ == "__main__":
    main()

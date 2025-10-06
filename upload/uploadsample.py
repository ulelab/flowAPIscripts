#!/usr/bin/env python3
"""
Flow upload_sample — multi-sample-type type_specific_metadata + post-upload GraphQL update.

- REST /upload/sample:
    * top-level: project, organism, sample_type (read from TSV 'Type')
    * type_specific_metadata (JSON): only ALLOWED keys for the inferred sample type
- Post-upload (GraphQL):
    * updateSample(...) via flat arguments on Mutation, selecting fields under updateSample { sample { ... } }
    * We DO NOT use purificationTargetId. We pass purificationTarget as a String (from TSV).
"""

import argparse
import csv
import getpass
import json
import os
import sys
import traceback
from typing import List, Dict, Any

import flowbio
import inspect  # optional; safe to keep


PROJECT_ID = 835660727529078273
INTERNAL_RETRIES = 5
DEFAULT_CHUNK_SIZE = 1_000_000

# ---- Per-sample-type whitelist of type_specific_metadata keys ----
# Sample type is read from TSV column 'Type' and used to filter keys below.
ALLOWED_TSM: Dict[str, set] = {
    "CLIP": {
        # barcodes/adapters/primers/UMIs that are typically safe for CLIP
        "five_prime_barcode_sequence",
        "three_prime_barcode_sequence",
        "three_prime_adapter_name",
        "three_prime_adapter_sequence",
        "rt_primer",
        "read1_primer",
        "read2_primer",
        "umi_barcode_sequence",
        "umi_separator",
        # (we intentionally omit experimental_method and purification_target for CLIP)
    },
    # RNA-Seq: strandedness and RNA selection are type-specific metadata
    "RNA-SEQ": {
        "strandedness",
        "rna_selection_method",
        # Some datasets may also include adapter/barcode fields; allow conservatively if present
        "three_prime_adapter_name",
        "three_prime_adapter_sequence",
        "read1_primer",
        "read2_primer",
    },
    # Ribo-Seq: ribosome and preparation-specific metadata
    "RIBO-SEQ": {
        "ribosome_type",
        "size_selection",
        "separation_method",
        "stabilization_method",
        # allow common library prep fields if provided
        "three_prime_adapter_name",
        "three_prime_adapter_sequence",
        "read1_primer",
        "read2_primer",
    },
}

# ---- Controlled vocabulary normalization -------------------------------------
def _norm_choice(value: str, mapping: Dict[str, str]) -> str:
    v = (value or "").strip().lower()
    if not v:
        return ""
    if v in mapping:
        return mapping[v]
    compact = "".join(ch for ch in v if ch.isalnum())
    for k, out in mapping.items():
        kcompact = "".join(ch for ch in k if ch.isalnum())
        if compact == kcompact:
            return out
    return value.strip()

RNA_STRANDEDNESS_MAP: Dict[str, str] = {
    "reverse": "reverse",
    "forward": "forward",
    "unstranded": "unstranded",
    "auto": "auto",
    # aliases often used in libraries/TSVs
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

# --- GraphQL: schema introspection (helpful once) -----------------------------
INTROSPECT_MUTATIONS = """
query Introspect {
  __schema {
    mutationType {
      fields {
        name
        args { name type { kind name ofType { kind name } } }
      }
    }
  }
}
"""

def print_sample_mutations(client):
    """Print mutation names/args that include 'sample' to help identify shapes."""
    try:
        resp = client.execute(INTROSPECT_MUTATIONS)
        fields = resp["data"]["__schema"]["mutationType"]["fields"]
        sample_fields = [f for f in fields if "sample" in f["name"].lower()]
        print("[DEBUG] Available Sample mutations & args:")
        for f in sample_fields:
            arg_list = []
            for a in f.get("args", []):
                t = a["type"]
                tname = t["name"] or (t["ofType"]["name"] if t.get("ofType") else None)
                arg_list.append(f"{a['name']}: {tname}")
            print(f"  - {f['name']}({', '.join(arg_list)})")
    except Exception as e:
        print(f"[DEBUG] Introspection failed (non-fatal): {e}")

# --- GraphQL update: ---
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
    sample {
      id
      name
      private
      sampleType
      organism { id name }
      source { id name }
      project { id name }
      purificationTarget { id name }   # ← fix: object with subfields
      purificationTargetText
      purificationAgent
      experimentalMethod
      condition
      sequencer
      scientist
      pi
      organisation
      fivePrimeBarcodeSequence
      threePrimeBarcodeSequence
      threePrimeAdapterName
      threePrimeAdapterSequence
      read1Primer
      read2Primer
      rtPrimer
      umiBarcodeSequence
      umiSeparator
      strandedness
      rnaSelectionMethod
      geo
      ena
      pubmed
      comments
    }
  }
}
"""


def build_update_vars(sample_id: str, meta: Dict[str, Any]) -> Dict[str, Any]:
    """
    Map TSV-derived metadata to Flow.bio's accepted camelCase fields.
    Only include non-empty values. NOTE: purificationTarget is a String (no IDs here).
    """
    alias_map = {
        "name": ["name", "Sample Name"],
        "organism": ["organism", "Organism"],
        "source": ["source", "Cell or Tissue", "Source"],
        "experimentalMethod": ["experimentalMethod", "Experimental Method", "experimental_method", "method"],
        "purificationAgent": ["purificationAgent", "Purification Agent", "purification_agent"],
        "purificationTarget": ["purificationTarget", "Protein (Purification Target)", "purification_target"],
        "purificationTargetText": ["purificationTargetText", "Purification Target Annotation", "purification_target_text"],
        "condition": ["condition", "Condition"],
        "sequencer": ["sequencer", "Sequencer", "platform"],
        "scientist": ["scientist", "Scientist"],
        "pi": ["pi", "PI"],
        "organisation": ["organisation", "Organisation", "organization", "Institute"],
        "fivePrimeBarcodeSequence": ["fivePrimeBarcodeSequence", "5' Barcode Sequence", "five_prime_barcode_sequence"],
        "threePrimeBarcodeSequence": ["threePrimeBarcodeSequence", "3' Barcode Sequence", "three_prime_barcode_sequence"],
        "threePrimeAdapterName": ["threePrimeAdapterName", "3' Adapter Name", "three_prime_adapter_name"],
        "threePrimeAdapterSequence": ["threePrimeAdapterSequence", "3' Adapter Sequence", "three_prime_adapter_sequence"],
        "read1Primer": ["read1Primer", "Read 1 Primer", "read1_primer"],
        "read2Primer": ["read2Primer", "Read 2 Primer", "read2_primer"],
        "rtPrimer": ["rtPrimer", "RT Primer", "rt_primer"],
        "umiBarcodeSequence": ["umiBarcodeSequence", "UMI Barcode Sequence", "umi_barcode_sequence"],
        "umiSeparator": ["umiSeparator", "UMI Separator", "umi_separator"],
        "strandedness": ["strandedness", "Strandedness"],
        "rnaSelectionMethod": ["rnaSelectionMethod", "RNA Selection Method", "rna_selection_method"],
        "geo": ["geo", "GEO ID", "GEO"],
        "ena": ["ena", "ENA ID", "ENA"],
        "pubmed": ["pubmed", "PubMed ID", "PMID"],
        "comments": ["comments", "Comments", "Notes"],
        "project": ["project", "Project"],  # usually avoid changing
        "private": ["private"],
    }

    def first_present(keys):
        for k in keys:
            if k in meta and str(meta[k]).strip():
                return str(meta[k]).strip()
        return None

    vars_out = {"id": str(sample_id)}
    for canon, candidates in alias_map.items():
        val = first_present(candidates)
        if val is not None:
            if canon == "private":
                val = True if val.lower() in {"1", "true", "yes"} else False
            vars_out[canon] = val

    # Avoid changing project here unless you *intend* to
    vars_out.pop("project", None)

    # DO NOT send empty strings
    return {k: v for k, v in vars_out.items() if v not in ("", None)}

def update_sample_metadata(client, sample_id: str, meta_from_tsv: Dict[str, Any]):
    """
    Call the fixed updateSample mutation with variables built from TSV.
    """
    variables = build_update_vars(sample_id, meta_from_tsv)
    if len(variables) == 1:  # only 'id'
        print(f"[DEBUG] No post-upload fields to update for sample {sample_id}.")
        return None
    try:
        resp = client.execute(UPDATE_SAMPLE_MUTATION, variables=variables)
        print("[DEBUG] GraphQL updateSample response:", resp)
        return resp["data"]["updateSample"]["sample"]
    except Exception as e:
        print(f"[ERROR] GraphQL updateSample failed for {sample_id}: {e}", file=sys.stderr)
        return None

# ---- Monkey-patch upload_sample to log first POST response -------------------
def _patched_upload_sample(self, name, path1, path2=None, chunk_size=1_000_000, progress=False,
                           metadata=None, use_base64=False, retries=0):
    import io, os, math, base64, requests
    from pathlib import Path
    from tqdm import tqdm

    class _Temp(io.BytesIO):
        def __init__(self, *args, name="", **kwargs):
            self.name = name
            super().__init__(*args, **kwargs)

    reads = [path1, path2] if path2 else [path1]
    data_id, sample_id, previous_data = None, None, []
    first_logged = False

    for path in reads:
        size = os.path.getsize(path)
        chunks = math.ceil(size / chunk_size)
        chunk_nums = tqdm(range(chunks)) if progress else range(chunks)
        for chunk_num in chunk_nums:
            filename = Path(path).name
            if progress:
                chunk_nums.set_description(f"Uploading {filename}")
            with open(path, "rb") as f:
                f.seek(chunk_num * chunk_size)
                data = f.read(chunk_size)
                if use_base64:
                    data = base64.b64encode(data)
                blob = _Temp(data, name=filename)

            is_last_data = (chunk_num == chunks - 1)
            is_last_sample = is_last_data and path == reads[-1]
            form = {
                "filename": filename,
                "is_last_sample": is_last_sample,
                "is_last": is_last_data,
                "expected_file_size": chunk_num * chunk_size,
                "sample_name": name,
                "data": data_id,
                "previous_data": previous_data,
                **(metadata or {}),
            }
            resp = requests.post(
                self.url.replace("/graphql", "/upload/sample"),
                data=form,
                files={"blob": blob},
                headers={"Authorization": self.headers["Authorization"]},
            )

            if (not first_logged) or (resp.status_code != 200):
                first_logged = True
                print(f"[DEBUG] upload_sample POST status={resp.status_code}")
                try:
                    print("[DEBUG] upload_sample JSON:", resp.json())
                except Exception:
                    print("[DEBUG] upload_sample TEXT:", resp.text[:2000])

            resp.raise_for_status()
            j = resp.json()
            if "data_id" not in j or "sample_id" not in j:
                raise RuntimeError(f"upload_sample: missing keys in response JSON: {j}")

            data_id = j["data_id"]
            sample_id = j["sample_id"]

            if is_last_data:
                previous_data.append(data_id)
                data_id = None

    return self.sample(sample_id)

flowbio.upload.UploadClient.upload_sample = _patched_upload_sample

# ---- Helpers -----------------------------------------------------------------
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

def add_if(d: Dict[str, Any], k: str, v: str):
    v = (v or "").strip()
    if v:
        d[k] = v

# ---- Main --------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Flow upload from TSV (multi-sample-type + GraphQL post-update)")
    ap.add_argument("tsv", help="Tab-delimited TSV")
    ap.add_argument("--rows", required=True, help="1-based rows, e.g. 1-10,15,22-24")
    ap.add_argument("--base-dir", default=".", help="Base dir for relative File paths")
    ap.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE)
    args = ap.parse_args()

    # TSV
    try:
        with open(args.tsv, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
    except Exception:
        print("Failed to read TSV", file=sys.stderr)
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

    # Discover the GraphQL execute method (usually client.execute)
    if not hasattr(client, "execute"):
        # Fall back: try to locate a callable attribute that looks like an executor
        candidates = [name for name, obj in inspect.getmembers(client) if callable(obj) and name in {"execute", "graphql", "gql"}]
        if candidates:
            # Alias it for our use below
            setattr(client, "execute", getattr(client, candidates[0]))
            print(f"[DEBUG] Using client.{candidates[0]} as GraphQL executor.")
        else:
            print("[ERROR] Could not find a GraphQL executor (expected client.execute).", file=sys.stderr)
            sys.exit(2)

    # Print available sample mutations once (diagnostic)
    print_sample_mutations(client)

    for idx in sel:
        row = rows[idx - 1]
        sample_name = (row.get("Sample Name") or "").strip()
        file_col    = (row.get("File") or "").strip()
        organism    = (row.get("Organism") or "").strip()
        sample_type = (row.get("Type") or "").strip()

        if not sample_name:
            print(f"[WARN] Row {idx}: missing 'Sample Name'; skipping.", file=sys.stderr)
            continue
        if not file_col:
            print(f"[WARN] Row {idx} ({sample_name}): missing 'File'; skipping.", file=sys.stderr)
            continue
        if not organism:
            print(f"[WARN] Row {idx} ({sample_name}): missing 'Organism'; skipping.", file=sys.stderr)
            continue

        path1 = norm_path(file_col, args.base_dir)
        if not os.path.exists(path1):
            print(f"[WARN] Row {idx} ({sample_name}): file not found: {path1}; skipping.", file=sys.stderr)
            continue

        # ---------- Build type_specific_metadata with whitelist ----------
        tsm: Dict[str, Any] = {
            "five_prime_barcode_sequence": (row.get("5' Barcode Sequence") or "").strip(),
            "three_prime_barcode_sequence": (row.get("3' Barcode Sequence") or "").strip(),
            "three_prime_adapter_name": (row.get("3' Adapter Name") or "").strip(),
            "three_prime_adapter_sequence": (row.get("3' Adapter Sequence") or "").strip(),
            "rt_primer": (row.get("RT Primer") or "").strip(),
            "read1_primer": (row.get("Read 1 Primer") or "").strip(),
            "read2_primer": (row.get("Read 2 Primer") or "").strip(),
            "umi_barcode_sequence": (row.get("UMI Barcode Sequence") or "").strip(),
            "umi_separator": (row.get("UMI Separator") or "").strip(),
            # Normalize RNA-Seq controlled vocabularies
            "strandedness": _norm_choice(row.get("Strandedness") or "", RNA_STRANDEDNESS_MAP),
            "rna_selection_method": _norm_choice(row.get("RNA Selection Method") or "", RNA_SELECTION_MAP),
            # Ribo-Seq: standardize size selection key name
            "size_selection": _norm_choice((row.get("Size selection") or row.get("Size Selection") or ""), RIBO_SIZE_SELECTION_MAP),
            "ribosome_type": _norm_choice((row.get("Ribosome Type") or ""), RIBO_TYPE_MAP),
            "separation_method": _norm_choice((row.get("Separation Method") or ""), RIBO_SEPARATION_METHOD_MAP),
            "stabilization_method": (row.get("Ribosome stabilization method") or "").strip(),
            # DO NOT include experimental_method / purification_target for CLIP (causes 400)
        }
        # Drop blanks
        tsm = {k: v for k, v in tsm.items() if v}

        # Filter by whitelist for this sample type (case-insensitive match)
        st_key = (sample_type or "").strip().upper()
        allowed = ALLOWED_TSM.get(st_key, set())
        if allowed:
            tsm = {k: v for k, v in tsm.items() if k in allowed}
        else:
            # If sample type is unknown, keep only very safe library-prep keys
            safe_keys = {
                "three_prime_adapter_name",
                "three_prime_adapter_sequence",
                "read1_primer",
                "read2_primer",
            }
            tsm = {k: v for k, v in tsm.items() if k in safe_keys}

        # Soft validation/warnings for required or expected fields per sample type
        def warn(msg: str):
            print(f"[WARN] Row {idx} ({sample_name}) [{sample_type or 'UNKNOWN'}]: {msg}", file=sys.stderr)

        if st_key == "RNA-SEQ":
            if not ((row.get("Strandedness") or "").strip()):
                warn("Missing Strandedness for RNA-Seq.")
            if not ((row.get("RNA Selection Method") or "").strip()):
                warn("Missing RNA Selection Method for RNA-Seq.")
        elif st_key == "RIBO-SEQ":
            if not ((row.get("Ribosome Type") or "").strip()):
                warn("Missing Ribosome Type for Ribo-Seq.")
            if not ((row.get("Size selection") or row.get("Size Selection") or "").strip()):
                warn("Missing Size selection for Ribo-Seq.")
            if not ((row.get("Separation Method") or "").strip()):
                warn("Missing Separation Method for Ribo-Seq.")
            if not ((row.get("Ribosome stabilization method") or "").strip()):
                warn("Missing Ribosome stabilization method for Ribo-Seq.")

        # ---------- Top-level metadata for REST upload ----------
        metadata: Dict[str, Any] = {
            "project": PROJECT_ID,
            "organism": organism,
            "sample_type": sample_type,  # as per docs
        }
        if tsm:
            metadata["type_specific_metadata"] = json.dumps(tsm)

        # Generic fields that have persisted for you
        add_if(metadata, "pi", row.get("PI", ""))
        add_if(metadata, "scientist", row.get("Scientist", ""))
        add_if(metadata, "organisation", row.get("Organisation", ""))
        add_if(metadata, "purificationAgent", row.get("Purification Agent", ""))
        add_if(metadata, "condition", row.get("Condition", ""))
        add_if(metadata, "sequencer", row.get("Sequencer", ""))
        add_if(metadata, "comments", row.get("Comments", ""))
        add_if(metadata, "source", row.get("Cell or Tissue", ""))
        add_if(metadata, "geo", row.get("GEO ID", ""))
        add_if(metadata, "ena", row.get("ENA ID", ""))
        add_if(metadata, "pubmed", row.get("PubMed ID", ""))

        # Strip empties
        metadata = {k: v for k, v in metadata.items() if v not in ("", None)}

        print(f"[DEBUG] Row {idx} '{sample_name}' outgoing metadata:", metadata)
        if "type_specific_metadata" in metadata:
            try:
                print("[DEBUG] type_specific_metadata (parsed):", json.loads(metadata["type_specific_metadata"]))
            except Exception:
                pass

        try:
            s = client.upload_sample(
                sample_name,
                path1,
                progress=True,
                retries=INTERNAL_RETRIES,
                chunk_size=args.chunk_size,
                metadata=metadata,
            )

            payload = s if isinstance(s, dict) else vars(s) if hasattr(s, "__dict__") else repr(s)
            sid = payload.get("id") if isinstance(payload, dict) else getattr(s, "id", None)

            print(f"[OK] Uploaded row {idx}: '{sample_name}' (id={sid})")
            print(f"[DEBUG] Row {idx} API response dump:", payload)

            # -------- Post-upload GraphQL update (no IDs required) --------
            # Build a TSV-shaped dict and hand to build_update_vars:
            tsv_meta = dict(row)  # shallow copy of CSV row
            # Also ensure we have the key we want for purificationTarget text
            if "Protein (Purification Target)" not in tsv_meta and "Purification Target" in tsv_meta:
                tsv_meta["Protein (Purification Target)"] = (row.get("Purification Target") or "").strip()

            if sid:
                update_sample_metadata(client, sid, tsv_meta)

        except Exception as e:
            print(f"[ERROR] Upload failed for row {idx} '{sample_name}': {e}", file=sys.stderr)
            traceback.print_exc()

    print("[DONE]")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# This helper is used to repeatedly convert Flow project JSON exports into a clean
# data_id list for bulk FASTQ download jobs (e.g. Slurm), so the same workflow can
# be rerun for different projects without manual copy/paste.

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List, Set


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract unique uploaded data IDs from flow_sample_data_pull JSON output "
            "and write one ID per line for bulk download scripts."
        )
    )
    parser.add_argument(
        "-i",
        "--input-json",
        required=True,
        help="Path to project_*_uploaded_data.json",
    )
    parser.add_argument(
        "-o",
        "--output-txt",
        required=True,
        help="Path to output data_id.txt (one ID per line)",
    )
    parser.add_argument(
        "--sort",
        action="store_true",
        help="Sort IDs before writing output",
    )
    return parser.parse_args()


def load_project_records(path: Path) -> List[Dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, list):
        raise ValueError(f"Expected top-level list in JSON: {path}")
    return [x for x in payload if isinstance(x, dict)]


def collect_uploaded_data_ids(records: List[Dict[str, Any]]) -> List[str]:
    unique: Set[str] = set()
    ordered: List[str] = []
    for record in records:
        uploaded = record.get("uploaded_data", [])
        if not isinstance(uploaded, list):
            continue
        for entry in uploaded:
            if not isinstance(entry, dict):
                continue
            data_id = entry.get("id")
            if not data_id:
                continue
            sid = str(data_id).strip()
            if sid and sid not in unique:
                unique.add(sid)
                ordered.append(sid)
    return ordered


def main() -> None:
    args = parse_args()
    input_path = Path(args.input_json)
    output_path = Path(args.output_txt)

    records = load_project_records(input_path)
    data_ids = collect_uploaded_data_ids(records)
    if args.sort:
        data_ids = sorted(data_ids)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(data_ids) + ("\n" if data_ids else ""), encoding="utf-8")

    print(
        f"Extracted {len(data_ids)} uploaded data IDs from {input_path} -> {output_path}"
    )


if __name__ == "__main__":
    main()

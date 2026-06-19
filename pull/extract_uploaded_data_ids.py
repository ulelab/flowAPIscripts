#!/usr/bin/env python3
# This helper is used to repeatedly convert Flow project JSON exports into a clean
# data_id list for bulk FASTQ download jobs (e.g. Slurm), so the same workflow can
# be rerun for different projects without manual copy/paste.

import argparse
import json
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Pattern, Set


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
    parser.add_argument(
        "--source-field",
        choices=("uploaded_data", "data"),
        default="uploaded_data",
        help="Which per-sample list in the JSON to read IDs from (default: uploaded_data)",
    )
    parser.add_argument(
        "--filename-regex",
        default=None,
        help="Optional extra filename filter when extracting IDs",
    )
    return parser.parse_args()


def load_project_records(path: Path) -> List[Dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, list):
        raise ValueError(f"Expected top-level list in JSON: {path}")
    return [x for x in payload if isinstance(x, dict)]


def _matches_filename(entry: Dict[str, Any], pattern: Optional[Pattern[str]]) -> bool:
    if pattern is None:
        return True
    filename = entry.get("filename") or ""
    return bool(pattern.search(str(filename)))


def collect_data_ids(
    records: List[Dict[str, Any]],
    source_field: str,
    filename_pattern: Optional[Pattern[str]] = None,
) -> List[str]:
    unique: Set[str] = set()
    ordered: List[str] = []
    for record in records:
        entries = record.get(source_field, [])
        if not isinstance(entries, list):
            continue
        for entry in entries:
            if not isinstance(entry, dict):
                continue
            if not _matches_filename(entry, filename_pattern):
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

    filename_pattern = re.compile(args.filename_regex) if args.filename_regex else None
    records = load_project_records(input_path)
    data_ids = collect_data_ids(
        records,
        source_field=args.source_field,
        filename_pattern=filename_pattern,
    )
    if args.sort:
        data_ids = sorted(data_ids)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(data_ids) + ("\n" if data_ids else ""), encoding="utf-8")

    print(
        f"Extracted {len(data_ids)} data IDs from {input_path} "
        f"(field={args.source_field}) -> {output_path}"
    )


if __name__ == "__main__":
    main()

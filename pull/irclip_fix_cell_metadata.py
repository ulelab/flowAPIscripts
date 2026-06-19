#!/usr/bin/env python3
"""
Build irCLIP cell-line fix CSV using original sample names as authoritative.

Authors confirmed GEO characteristics are wrong; sample names (from audit pull)
are the source of truth. Also applies A431 -> HEK293T renames where the
original upload name used HEK293T but Flow was later changed to A431.
"""

from __future__ import annotations

import argparse
import csv
import re
from copy import deepcopy
from pathlib import Path

# Original LOW-confidence HNRNPC rows: authors / user want HEK293T not A431.
FORCE_HEK293T_SAMPLE_IDS = {
    "118950627770000171",
    "421233802465164873",
    "822557063120492700",
}


def cell_from_name(name: str) -> str:
    m = re.search(r"_Hs_([^_]+)_", name)
    return m.group(1) if m else ""


def authoritative_name(audit_name: str, sample_id: str) -> str:
    if sample_id in FORCE_HEK293T_SAMPLE_IDS:
        return audit_name.replace("_Hs_A431_", "_Hs_HEK293T_")
    return audit_name


def apply_fixes(row: dict[str, str], audit_name: str) -> tuple[dict[str, str], list[str]]:
    out = deepcopy(row)
    notes: list[str] = []
    sid = out["sample_id"]
    target_name = authoritative_name(audit_name, sid)
    target_source = cell_from_name(target_name)

    if out["sample_name"] != target_name:
        notes.append(f"rename:{out['sample_name']}->{target_name}")
        out["sample_name"] = target_name

    if target_source and out.get("source", "") != target_source:
        notes.append(f"source:{out.get('source', '')}->{target_source}")
        out["source"] = target_source

    return out, notes


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline", required=True, help="Current Flow pull CSV")
    ap.add_argument("--audit", required=True, help="Original sample_name_audit_full.csv")
    ap.add_argument("--output-updated", required=True)
    ap.add_argument("--output-manifest", required=True)
    args = ap.parse_args()

    audit_by_id = {r["sample_id"]: r["sample_name"] for r in csv.DictReader(open(args.audit, newline="", encoding="utf-8"))}
    rows = list(csv.DictReader(open(args.baseline, newline="", encoding="utf-8")))
    changed_rows: list[dict[str, str]] = []
    manifest: list[dict[str, str]] = []

    for row in rows:
        sid = row["sample_id"]
        audit_name = audit_by_id.get(sid, row["sample_name"])
        fixed, notes = apply_fixes(row, audit_name)
        if notes:
            changed_rows.append(fixed)
            manifest.append(
                {
                    "sample_id": sid,
                    "project_id": row.get("project_id", ""),
                    "authoritative_sample_name": authoritative_name(audit_name, sid),
                    "old_sample_name": row["sample_name"],
                    "new_sample_name": fixed["sample_name"],
                    "old_source": row.get("source", ""),
                    "new_source": fixed.get("source", ""),
                    "cell_from_name": cell_from_name(fixed["sample_name"]),
                    "changes": "; ".join(notes),
                    "flow_sample_url": f"https://app.flow.bio/samples/{sid}",
                }
            )

    if changed_rows:
        with open(args.output_updated, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(changed_rows[0].keys()))
            w.writeheader()
            w.writerows(changed_rows)

    if manifest:
        with open(args.output_manifest, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(manifest[0].keys()))
            w.writeheader()
            w.writerows(manifest)

    print(f"Changed samples: {len(changed_rows)}")
    for m in manifest:
        print(f"  {m['old_sample_name']}")
        if m["old_sample_name"] != m["new_sample_name"]:
            print(f"    name -> {m['new_sample_name']}")
        if m["old_source"] != m["new_source"]:
            print(f"    source {m['old_source']} -> {m['new_source']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

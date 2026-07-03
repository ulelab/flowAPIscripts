#!/usr/bin/env python3
"""Merge agent-confirmed metadata proposals into an updated pull CSV for push v2."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def _norm(v: object) -> str:
    if v is None:
        return ""
    s = str(v).strip()
    return "" if s.lower() in ("nan", "none") else s


def load_csv(path: Path) -> dict[str, dict[str, str]]:
    by_id: dict[str, dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            sid = _norm(row.get("sample_id"))
            if sid:
                by_id[sid] = {k: _norm(v) for k, v in row.items()}
    return by_id


def load_proposals(path: Path) -> tuple[dict[str, str], dict[str, str], list[str]]:
    """Return target->value, sample_id->value overrides, evidence lines."""
    by_target: dict[str, str] = {}
    by_sample: dict[str, str] = {}
    evidence: list[str] = []
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames and "purification_target" not in (reader.fieldnames or []):
            f.seek(0)
            reader = csv.DictReader(f)
        for row in reader:
            status = _norm(row.get("status", "pending_confirmation"))
            if status != "confirmed":
                continue
            val = _norm(row.get("proposed_purification_agent") or row.get("proposed_value"))
            sid = _norm(row.get("sample_id"))
            target = _norm(row.get("purification_target"))
            quote = _norm(row.get("evidence_quote"))
            if sid:
                by_sample[sid] = val
            elif target:
                by_target[target] = val
            if quote and target:
                evidence.append(f"- **{target}**: {val!r} — {quote}")
    return by_target, by_sample, evidence


def write_confirm_md(
    path: Path,
    *,
    field: str,
    changes: list[tuple[str, str, str, str]],
    evidence: list[str],
) -> None:
    lines = [
        f"# Metadata update review — `{field}`",
        "",
        "Human must set `status: confirmed` on proposal rows before running this script.",
        "",
        "## Changes",
        "",
        "| sample_id | sample_name | old | new |",
        "|-----------|-------------|-----|-----|",
    ]
    for sid, name, old, new in changes:
        lines.append(f"| {sid} | {name} | {old!r} | {new!r} |")
    if evidence:
        lines.extend(["", "## Evidence (from proposals)", ""] + evidence)
    lines.extend(
        [
            "",
            "## Push",
            "",
            "```bash",
            "python3 flowAPIscripts/pull/flow_public_samples_push_metadata_v2.py \\",
            "  --dry-run --baseline samples_baseline.csv --updated samples_updated.csv",
            "```",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--baseline", type=Path, required=True)
    ap.add_argument("--proposals", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--field", default="purification_agent")
    ap.add_argument(
        "--confirm-md",
        type=Path,
        default=None,
        help="Write human-readable change summary (default: output dir CONFIRM_METADATA_UPDATES.md)",
    )
    args = ap.parse_args()

    baseline = load_csv(args.baseline)
    by_target, by_sample, evidence = load_proposals(args.proposals)
    if not by_target and not by_sample:
        raise SystemExit("No confirmed proposals found (status must be 'confirmed').")

    updated = {sid: dict(row) for sid, row in baseline.items()}
    changes: list[tuple[str, str, str, str]] = []

    for sid, row in updated.items():
        target = _norm(row.get("purification_target"))
        if target == "SMInput":
            continue
        new_val = by_sample.get(sid)
        if new_val is None and target in by_target:
            new_val = by_target[target]
        if new_val is None:
            continue
        old_val = _norm(row.get(args.field))
        if old_val == new_val:
            continue
        row[args.field] = new_val
        changes.append((sid, _norm(row.get("sample_name")), old_val, new_val))

    if not changes:
        raise SystemExit("No field changes after applying confirmed proposals.")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = sorted({k for row in updated.values() for k in row})
    with args.output.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for sid in sorted(updated):
            w.writerow(updated[sid])

    confirm = args.confirm_md or (args.output.parent / "CONFIRM_METADATA_UPDATES.md")
    write_confirm_md(confirm, field=args.field, changes=changes, evidence=evidence)
    print(f"Wrote {args.output} ({len(changes)} sample(s) changed)")
    print(f"Review {confirm}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

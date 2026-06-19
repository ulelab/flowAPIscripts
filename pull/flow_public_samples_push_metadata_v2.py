#!/usr/bin/env python3
"""
Diff updated public-samples CSV against a baseline pull, then push metadata to Flow.bio.

v2 adds REST /edit for attribute annotations on nested metadata objects:
  - purification_target + purification_target__annotation
  - source + source__annotation (annotation is the sub-field on source)

Transport:
  - GraphQL updateSample: scalar metadata fields (not *_ _annotation)
  - REST POST app.flow.bio/api/samples/{id}/edit: *_ _annotation keys (flat; matches GET)

Credentials: FLOWBIO_USERNAME / FLOWBIO_PASSWORD or --username / --password.

Examples:
  python3 flow_public_samples_push_metadata_v2.py --dry-run \\
    --baseline flow_public_samples_pull_v7.csv \\
    --updated flow_public_samples_bulk_push_w7_colleague_updates.csv

  python3 flow_public_samples_push_metadata_v2.py --yes --allow-clear \\
    --baseline flow_public_samples_pull_v7.csv \\
    --updated flow_public_samples_bulk_push_w7_colleague_updates.csv \\
    --username USER --password PASS
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import sys
import warnings
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# Compatibility shim for Python versions without warnings.deprecated (flowbio import).
if not hasattr(warnings, "deprecated"):

    def _deprecated(_msg: str):  # type: ignore[no-redef]
        def _decorator(func):
            return func

        return _decorator

    warnings.deprecated = _deprecated  # type: ignore[attr-defined]

import requests
from flowbio import Client

# Reuse REST helpers from pull script (same directory).
_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from flow_public_samples_pull_v3 import (  # noqa: E402
    API_BASE,
    fetch_sample_detail,
    flatten_sample_detail,
    rest_login,
)

DEFAULT_BASELINE = _SCRIPT_DIR / "flow_public_samples_pull_v7.csv"
DEFAULT_UPDATED = _SCRIPT_DIR / "flow_public_samples_bulk_push_w7_colleague_updates.csv"
REST_EDIT_API_BASE = "https://app.flow.bio/api"

# REST /edit: annotation CSV column -> parent value column on same metadata object
REST_ANNOTATION_PARENT: Dict[str, str] = {
    "purification_target__annotation": "purification_target",
    "source__annotation": "source",
}

# CSV column -> GraphQL variable name (do not use *Text GraphQL vars; they do not persist)
WHITELIST_GRAPHQL: Dict[str, str] = {
    "sample_name": "name",
    "condition": "condition",
    "comments": "comments",
    "experimental_method": "experimentalMethod",
    "purification_agent": "purificationAgent",
    "purification_target": "purificationTarget",
    "source": "source",
}

# CSV columns pushed via REST /edit (flat keys; same names as pull CSV)
WHITELIST_REST_EDIT: Tuple[str, ...] = tuple(REST_ANNOTATION_PARENT.keys())

GRAPHQL_VAR_TYPES: Dict[str, str] = {
    "name": "String",
    "condition": "String",
    "comments": "String",
    "experimentalMethod": "String",
    "purificationAgent": "String",
    "purificationTarget": "String",
    "source": "String",
}


def _norm(value: Any) -> str:
    if value is None:
        return ""
    s = str(value).strip()
    if s.lower() in ("nan", "none"):
        return ""
    return s


def load_csv_indexed(path: Path) -> Dict[str, Dict[str, str]]:
    """Load CSV keyed by sample_id; raises on duplicate sample_id."""
    if not path.is_file():
        raise FileNotFoundError(f"CSV not found: {path}")

    by_id: Dict[str, Dict[str, str]] = {}
    dupes: List[Tuple[str, int]] = []

    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row_num, raw in enumerate(reader, start=2):
            row = {k: _norm(v) for k, v in raw.items()}
            sid = row.get("sample_id", "")
            if not sid:
                logging.warning("Row %d: missing sample_id, skipped", row_num)
                continue
            if sid in by_id:
                dupes.append((sid, row_num))
                continue
            by_id[sid] = row

    if dupes:
        lines = ", ".join(f"{sid} (row {r})" for sid, r in dupes[:10])
        raise RuntimeError(f"Duplicate sample_id in {path}: {lines}")

    return by_id


@dataclass
class FieldChange:
    csv_column: str
    transport: str  # "graphql" | "rest"
    api_field: str
    old_value: str
    new_value: str


@dataclass
class PendingChange:
    sample_id: str
    sample_name_v3: str
    sample_name_v4: str
    v4_row: Dict[str, str] = field(default_factory=dict)
    field_changes: List[FieldChange] = field(default_factory=list)

    def graphql_payload(self) -> Dict[str, Any]:
        out: Dict[str, Any] = {"id": self.sample_id}
        for fc in self.field_changes:
            if fc.transport == "graphql":
                out[fc.api_field] = fc.new_value
        return out

    def rest_edit_body(self) -> Dict[str, str]:
        """Body for POST /samples/{id}/edit (flat keys)."""
        body: Dict[str, str] = {}
        for fc in self.field_changes:
            if fc.transport == "rest":
                body[fc.api_field] = fc.new_value
        # REST /edit expects parent object value alongside annotation (see GET metadata shape).
        for ann_col, parent_col in REST_ANNOTATION_PARENT.items():
            if ann_col in body or parent_col in body:
                parent_val = _norm(self.v4_row.get(parent_col, ""))
                if parent_val:
                    body[parent_col] = parent_val
        return body

    def payload(self) -> Dict[str, Any]:
        """Legacy combined view for audit logs."""
        out: Dict[str, Any] = {"id": self.sample_id}
        out.update(self.graphql_payload())
        out.update(self.rest_edit_body())
        return out

    def expected_values(self) -> Dict[str, str]:
        """CSV column -> expected value after push (from v4)."""
        base = {
            "sample_name": self.sample_name_v4,
            "condition": "",
            "comments": "",
            "purification_agent": "",
            "purification_target": "",
            "purification_target__annotation": "",
            "source": "",
            "source__annotation": "",
        }
        for fc in self.field_changes:
            base[fc.csv_column] = fc.new_value
        return base


def _value_changed(old_v: str, new_v: str, allow_clear: bool) -> bool:
    if old_v == new_v:
        return False
    if not new_v and not allow_clear:
        return False
    return True


def compute_pending_changes(
    baseline: Dict[str, Dict[str, str]],
    updated: Dict[str, Dict[str, str]],
    allow_clear: bool = False,
) -> Tuple[List[PendingChange], List[str]]:
    pending: List[PendingChange] = []
    warnings: List[str] = []

    for sid, v4_row in updated.items():
        v3_row = baseline.get(sid)
        if not v3_row:
            warnings.append(f"sample_id {sid}: in updated file but not in baseline, skipped")
            continue

        field_changes: List[FieldChange] = []

        for csv_col in WHITELIST_REST_EDIT:
            old_v = _norm(v3_row.get(csv_col, ""))
            new_v = _norm(v4_row.get(csv_col, ""))
            if not _value_changed(old_v, new_v, allow_clear):
                continue
            field_changes.append(
                FieldChange(
                    csv_column=csv_col,
                    transport="rest",
                    api_field=csv_col,
                    old_value=old_v,
                    new_value=new_v,
                )
            )

        rest_cols = {fc.csv_column for fc in field_changes}
        rest_parents = {REST_ANNOTATION_PARENT[c] for c in rest_cols if c in REST_ANNOTATION_PARENT}

        for csv_col, gql_field in WHITELIST_GRAPHQL.items():
            old_v = _norm(v3_row.get(csv_col, ""))
            new_v = _norm(v4_row.get(csv_col, ""))
            if not _value_changed(old_v, new_v, allow_clear):
                continue
            # When annotation changes, push parent via REST alongside annotation.
            if csv_col in rest_parents:
                field_changes.append(
                    FieldChange(
                        csv_column=csv_col,
                        transport="rest",
                        api_field=csv_col,
                        old_value=old_v,
                        new_value=new_v,
                    )
                )
                continue
            field_changes.append(
                FieldChange(
                    csv_column=csv_col,
                    transport="graphql",
                    api_field=gql_field,
                    old_value=old_v,
                    new_value=new_v,
                )
            )

        if field_changes:
            pending.append(
                PendingChange(
                    sample_id=sid,
                    sample_name_v3=_norm(v3_row.get("sample_name", "")),
                    sample_name_v4=_norm(v4_row.get("sample_name", "")),
                    v4_row=dict(v4_row),
                    field_changes=field_changes,
                )
            )

    pending.sort(key=lambda p: p.sample_id)
    return pending, warnings


def write_pending_csv(path: Path, pending: List[PendingChange]) -> None:
    cols = [
        "sample_id",
        "sample_name_v3",
        "sample_name_v4",
        "csv_column",
        "transport",
        "api_field",
        "old_value",
        "new_value",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for p in pending:
            for fc in p.field_changes:
                w.writerow(
                    {
                        "sample_id": p.sample_id,
                        "sample_name_v3": p.sample_name_v3,
                        "sample_name_v4": p.sample_name_v4,
                        "csv_column": fc.csv_column,
                        "transport": fc.transport,
                        "api_field": fc.api_field,
                        "old_value": fc.old_value,
                        "new_value": fc.new_value,
                    }
                )


def build_update_mutation(graphql_fields: List[str]) -> str:
    if not graphql_fields:
        raise ValueError("No GraphQL fields to update")

    var_decls = ["$id: ID!"]
    mutation_args = ["id: $id"]
    for gf in graphql_fields:
        gql_type = GRAPHQL_VAR_TYPES.get(gf, "String")
        var_decls.append(f"${gf}: {gql_type}")
        mutation_args.append(f"{gf}: ${gf}")

    var_block = ",\n      ".join(var_decls)
    args_block = ",\n        ".join(mutation_args)

    return f"""
    mutation UpdateSample(
      {var_block}
    ) {{
      updateSample(
        {args_block}
      ) {{
        sample {{
          id
          name
        }}
      }}
    }}
    """


def ensure_graphql_execute(client: Client) -> None:
    if hasattr(client, "execute") and callable(client.execute):
        return
    for candidate in ("execute", "graphql", "gql"):
        func = getattr(client, candidate, None)
        if callable(func):
            setattr(client, "execute", func)
            return
    raise RuntimeError("flowbio Client has no GraphQL execute method")


def rest_edit_sample(
    session: requests.Session,
    token: str,
    sample_id: str,
    body: Dict[str, str],
) -> Tuple[bool, Optional[str]]:
    if not body:
        return True, None
    headers = {
        "Authorization": f"Bearer {token}",
        "Content-Type": "application/json",
    }
    url = f"{REST_EDIT_API_BASE}/samples/{sample_id}/edit"
    try:
        r = session.post(url, json=body, headers=headers, timeout=90)
    except Exception as exc:  # noqa: BLE001
        return False, str(exc)
    if r.status_code >= 400:
        return False, f"HTTP {r.status_code}: {r.text[:800]}"
    return True, None


def push_graphql_update(client: Client, pending: PendingChange) -> Tuple[bool, Optional[str]]:
    gql_fields = [fc.api_field for fc in pending.field_changes if fc.transport == "graphql"]
    if not gql_fields:
        return True, None
    payload = pending.graphql_payload()
    mutation = build_update_mutation(gql_fields)

    try:
        result = client.execute(mutation, variables=payload)
    except Exception as exc:  # noqa: BLE001
        return False, str(exc)

    if result and "errors" in result:
        return False, json.dumps(result["errors"], indent=2)

    if result and result.get("data", {}).get("updateSample"):
        return True, None

    if result and "updateSample" in result:
        return True, None

    return False, json.dumps(result, indent=2) if result else "empty response"


def push_update(
    client: Client,
    session: requests.Session,
    token: str,
    pending: PendingChange,
) -> Tuple[bool, Optional[str]]:
    rest_body = pending.rest_edit_body()
    ok, err = rest_edit_sample(session, token, pending.sample_id, rest_body)
    if not ok:
        return False, f"REST edit failed: {err}"
    return push_graphql_update(client, pending)


def verify_sample(
    session: requests.Session,
    token: str,
    pending: PendingChange,
) -> Tuple[bool, List[str]]:
    detail = fetch_sample_detail(session, token, pending.sample_id)
    flat = flatten_sample_detail(detail)
    expected = pending.expected_values()
    lines: List[str] = []
    ok = True

    for fc in pending.field_changes:
        live = _norm(flat.get(fc.csv_column, ""))
        exp = _norm(expected.get(fc.csv_column, ""))
        if live == exp:
            lines.append(f"  {fc.csv_column}: OK ({live!r})")
        else:
            ok = False
            lines.append(
                f"  {fc.csv_column}: MISMATCH expected={exp!r} live={live!r}"
            )

    return ok, lines


def print_change_summary(pending: PendingChange) -> None:
    print(f"\nsample_id: {pending.sample_id}")
    if pending.sample_name_v3 != pending.sample_name_v4:
        print(f"sample_name (v3 -> v4): {pending.sample_name_v3!r} -> {pending.sample_name_v4!r}")
    for fc in pending.field_changes:
        via = "REST /edit" if fc.transport == "rest" else "GraphQL"
        print(f"  [{via}] {fc.csv_column}: {fc.old_value!r} -> {fc.new_value!r}")
    rest_body = pending.rest_edit_body()
    if rest_body:
        print("REST POST body (app.flow.bio/api/samples/.../edit):")
        print(json.dumps(rest_body, indent=2))
    gql = pending.graphql_payload()
    if len(gql) > 1:
        print("GraphQL variables:")
        print(json.dumps(gql, indent=2))


def run_single_sample_test(
    session: requests.Session,
    token: str,
    sample_id: str,
    annotation: str,
    parent_value: str = "",
    field: str = "purification_target__annotation",
) -> int:
    """POST one REST annotation edit and verify via GET."""
    if field not in REST_ANNOTATION_PARENT:
        print(f"Unsupported --test-field {field!r}; choose from {list(REST_ANNOTATION_PARENT)}")
        return 2
    parent_col = REST_ANNOTATION_PARENT[field]
    body: Dict[str, str] = {field: annotation}
    detail = fetch_sample_detail(session, token, sample_id)
    flat = flatten_sample_detail(detail)
    pv = parent_value or _norm(flat.get(parent_col, ""))
    if pv:
        body[parent_col] = pv
    print(f"Current live {parent_col}: {pv!r}")
    print(f"Current live {field}: {_norm(flat.get(field, ''))!r}")

    print(f"\nPOST {REST_EDIT_API_BASE}/samples/{sample_id}/edit")
    print(json.dumps(body, indent=2))
    ok, err = rest_edit_sample(session, token, sample_id, body)
    if not ok:
        print(f"FAILED: {err}")
        return 1
    print("POST OK")

    flat = flatten_sample_detail(fetch_sample_detail(session, token, sample_id))
    live = _norm(flat.get(field, ""))
    if live == annotation:
        print(f"VERIFY OK: {field}={live!r}")
        return 0
    print(f"VERIFY MISMATCH: expected={annotation!r} live={live!r}")
    return 1


def append_audit(audit_path: Path, record: Dict[str, Any]) -> None:
    with audit_path.open("a", encoding="utf-8") as f:
        f.write(json.dumps(record, ensure_ascii=False) + "\n")


def prompt_action(auto_yes: bool) -> str:
    if auto_yes:
        return "y"
    while True:
        ans = input("Push this update? [y]es / [n]o skip / [q]uit: ").strip().lower()
        if ans in ("y", "yes", "n", "no", "q", "quit"):
            return "y" if ans in ("y", "yes") else ("q" if ans in ("q", "quit") else "n")
        print("Invalid input; enter y, n, or q.")


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    ap.add_argument("--updated", type=Path, default=DEFAULT_UPDATED)
    ap.add_argument("--username", default=os.environ.get("FLOWBIO_USERNAME", ""))
    ap.add_argument("--password", default=os.environ.get("FLOWBIO_PASSWORD", ""))
    ap.add_argument("--dry-run", action="store_true", help="List pending changes only")
    ap.add_argument("--limit", type=int, default=0, help="Process at most N pending changes")
    ap.add_argument("--yes", action="store_true", help="Auto-approve all pushes (after dry-run review)")
    ap.add_argument(
        "--on-verify-fail",
        choices=("stop", "continue"),
        default="stop",
        help="If verification fails after push (default: stop)",
    )
    ap.add_argument(
        "--audit-dir",
        type=Path,
        default=_SCRIPT_DIR,
        help="Directory for audit jsonl and pending csv",
    )
    ap.add_argument(
        "--test-sample-id",
        default="",
        help="If set, POST one REST edit for this sample and verify (skip bulk diff)",
    )
    ap.add_argument(
        "--test-annotation",
        default="reCLIP_hnRNPC",
        help="With --test-sample-id, annotation value to POST",
    )
    ap.add_argument(
        "--test-field",
        default="purification_target__annotation",
        choices=tuple(REST_ANNOTATION_PARENT.keys()),
        help="With --test-sample-id, which annotation column to test",
    )
    ap.add_argument(
        "--test-parent-value",
        default="",
        help="With --test-sample-id, parent value (e.g. source or purification_target); default from live sample",
    )
    ap.add_argument(
        "--allow-clear",
        action="store_true",
        help="Push empty values when baseline had text (clears condition/annotations)",
    )
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")

    if not args.username or not args.password:
        logging.error("Set FLOWBIO_USERNAME/FLOWBIO_PASSWORD or use --username/--password")
        return 2

    session = requests.Session()
    try:
        token = rest_login(session, args.username, args.password)
    except Exception as exc:  # noqa: BLE001
        logging.error("REST login failed: %s", exc)
        return 2

    if args.test_sample_id:
        return run_single_sample_test(
            session,
            token,
            args.test_sample_id.strip(),
            args.test_annotation.strip(),
            parent_value=args.test_parent_value.strip(),
            field=args.test_field,
        )

    try:
        baseline = load_csv_indexed(args.baseline)
        updated = load_csv_indexed(args.updated)
    except (FileNotFoundError, RuntimeError) as exc:
        logging.error("%s", exc)
        return 2

    pending, warnings = compute_pending_changes(baseline, updated, allow_clear=args.allow_clear)
    for w in warnings:
        logging.warning("%s", w)

    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    pending_path = args.audit_dir / f"flow_public_samples_push_v2_pending_{ts}.csv"
    write_pending_csv(pending_path, pending)
    logging.info("Wrote %d pending field changes to %s", sum(len(p.field_changes) for p in pending), pending_path)
    logging.info("Samples with changes: %d", len(pending))

    if args.limit > 0:
        pending = pending[: args.limit]

    if args.dry_run:
        for p in pending:
            print_change_summary(p)
        print(f"\nDry run complete: {len(pending)} sample(s) with pushable changes.")
        return 0

    if not pending:
        logging.info("No pushable changes found.")
        return 0

    audit_path = args.audit_dir / f"flow_public_samples_push_v2_audit_{ts}.jsonl"

    client = Client()
    try:
        client.login(args.username, args.password)
        logging.info("Logged in to Flow.bio (GraphQL client)")
    except Exception as exc:  # noqa: BLE001
        logging.error("Login failed: %s", exc)
        return 2

    try:
        ensure_graphql_execute(client)
    except RuntimeError as exc:
        logging.error("%s", exc)
        return 2

    pushed = 0
    skipped = 0

    for p in pending:
        print_change_summary(p)
        action = prompt_action(args.yes)

        if action == "q":
            append_audit(
                audit_path,
                {
                    "timestamp": datetime.now(timezone.utc).isoformat(),
                    "sample_id": p.sample_id,
                    "action": "quit",
                    "fields_changed": [fc.csv_column for fc in p.field_changes],
                },
            )
            logging.info("Quit by user after %d pushed, %d skipped", pushed, skipped)
            break

        if action == "n":
            skipped += 1
            append_audit(
                audit_path,
                {
                    "timestamp": datetime.now(timezone.utc).isoformat(),
                    "sample_id": p.sample_id,
                    "action": "skipped",
                    "fields_changed": [fc.csv_column for fc in p.field_changes],
                    "payload": p.payload(),
                },
            )
            continue

        ok, err = push_update(client, session, token, p)
        if not ok:
            append_audit(
                audit_path,
                {
                    "timestamp": datetime.now(timezone.utc).isoformat(),
                    "sample_id": p.sample_id,
                    "action": "push_failed",
                    "rest_edit_body": p.rest_edit_body(),
                    "graphql_payload": p.graphql_payload(),
                    "error": err,
                },
            )
            logging.error("Push failed for %s: %s", p.sample_id, err)
            retry = input("Push failed. [r]etry / [n]skip / [q]uit: ").strip().lower()
            if retry in ("r", "retry"):
                ok, err = push_update(client, session, token, p)
                if not ok:
                    logging.error("Retry failed: %s", err)
                    if args.on_verify_fail == "stop":
                        return 1
                    skipped += 1
                    continue
            elif retry in ("q", "quit"):
                return 1
            else:
                skipped += 1
                continue

        verify_ok, verify_lines = verify_sample(session, token, p)
        print("Verification:")
        for line in verify_lines:
            print(line)

        append_audit(
            audit_path,
            {
                "timestamp": datetime.now(timezone.utc).isoformat(),
                "sample_id": p.sample_id,
                "action": "pushed",
                "fields_changed": [fc.csv_column for fc in p.field_changes],
                "rest_edit_body": p.rest_edit_body(),
                "graphql_payload": p.graphql_payload(),
                "verify_ok": verify_ok,
                "verify_lines": verify_lines,
            },
        )

        if not verify_ok:
            logging.error("Verification failed for sample_id %s", p.sample_id)
            if args.on_verify_fail == "stop":
                return 1
        else:
            pushed += 1
            logging.info("Pushed and verified sample_id %s", p.sample_id)

    logging.info("Done. pushed=%d skipped=%d audit=%s", pushed, skipped, audit_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())

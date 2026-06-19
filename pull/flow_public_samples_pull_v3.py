#!/usr/bin/env python3
"""
Pull full existing metadata for all public Flow.bio samples (private == False).

Uses the REST attribute shape: each metadata entry is typically
  { "value": ..., "annotation": ..., "attribute_name": ... }
Primary values are written under the metadata key (e.g. purification_target).
Secondary free-text slots are written as <key>__annotation (e.g. purification_target__annotation).

Requires authentication to list projects and samples. Public sample details can
also be fetched without auth, but discovery still needs /projects + /samples listing.

Credentials: FLOWBIO_USERNAME / FLOWBIO_PASSWORD or --username / --password.

Example:
  export FLOWBIO_USERNAME=...
  export FLOWBIO_PASSWORD=...
  python3 flow_public_samples_pull_v3.py --output-csv /home/mikej10/advbfx/projects/flow_public_samples_pull_v3.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

import requests

API_BASE = "https://api.flow.bio"
DEFAULT_WORKERS = 24

# Columns that always appear first if present
CORE_ORDER = [
    "project_id",
    "project_name",
    "sample_id",
    "sample_name",
    "private",
    "sample_type",
    "organism",
    "pubmed",
    "file_names",
]


def _raise(resp: requests.Response) -> None:
    try:
        resp.raise_for_status()
    except requests.HTTPError as e:
        raise requests.HTTPError(f"HTTP {resp.status_code}: {resp.text[:800]}") from e


def rest_login(session: requests.Session, username: str, password: str) -> str:
    r = session.post(
        f"{API_BASE}/login",
        json={"username": username, "password": password},
        timeout=90,
    )
    _raise(r)
    tok = r.json().get("token")
    if not tok:
        raise RuntimeError("Login failed: no token")
    return str(tok)


def parse_projects_payload(data: Any) -> List[Dict[str, Any]]:
    if isinstance(data, list):
        return data
    if isinstance(data, dict):
        for k in ("projects", "data", "items", "results"):
            v = data.get(k)
            if isinstance(v, list):
                return v
    raise RuntimeError(f"Unexpected /projects JSON shape: {type(data)}")


def fetch_all_projects(session: requests.Session, token: str, page_size: int = 100) -> List[Dict[str, Any]]:
    headers = {"Authorization": f"Bearer {token}"}
    page = 1
    out: List[Dict[str, Any]] = []
    while True:
        r = session.get(
            f"{API_BASE}/projects",
            params={"page": page, "count": page_size},
            headers=headers,
            timeout=90,
        )
        _raise(r)
        chunk = parse_projects_payload(r.json())
        out.extend(chunk)
        if len(chunk) < page_size:
            break
        page += 1
    return out


def fetch_project_samples_list(
    session: requests.Session, token: str, project_id: str, page_size: int = 100
) -> List[Dict[str, Any]]:
    headers = {"Authorization": f"Bearer {token}"}
    page = 1
    samples: List[Dict[str, Any]] = []
    while True:
        r = session.get(
            f"{API_BASE}/projects/{project_id}/samples",
            params={"page": page, "count": page_size},
            headers=headers,
            timeout=90,
        )
        _raise(r)
        chunk = r.json().get("samples") or []
        samples.extend(chunk)
        if len(chunk) < page_size:
            break
        page += 1
    return samples


def fetch_sample_detail(
    session: requests.Session, token: Optional[str], sample_id: str
) -> Dict[str, Any]:
    headers = {}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    r = session.get(f"{API_BASE}/samples/{sample_id}", headers=headers, timeout=90)
    _raise(r)
    return r.json()


def file_names_from_detail(detail: Dict[str, Any]) -> str:
    names: List[str] = []
    for fs in detail.get("filesets") or []:
        if not isinstance(fs, dict):
            continue
        for data in fs.get("data") or []:
            if isinstance(data, dict) and data.get("filename"):
                names.append(str(data["filename"]))
    return ",".join(names)


def flatten_sample_detail(detail: Dict[str, Any]) -> Dict[str, str]:
    """Turn GET /samples/{id} into a flat string row."""
    proj = detail.get("project") or {}
    if not isinstance(proj, dict):
        proj = {}

    row: Dict[str, str] = {
        "project_id": str(proj.get("id") or ""),
        "project_name": str(proj.get("name") or ""),
        "sample_id": str(detail.get("id") or ""),
        "sample_name": str(detail.get("name") or ""),
        "private": "" if detail.get("private") is None else str(detail.get("private")).lower(),
        "sample_type": str(detail.get("sample_type") or ""),
        "organism": str(detail.get("organism") or ""),
        "pubmed": str(detail.get("pubmed") or ""),
        "file_names": file_names_from_detail(detail),
    }

    md = detail.get("metadata") or {}
    if not isinstance(md, dict):
        return row

    for key, obj in md.items():
        if not isinstance(obj, dict):
            row[str(key)] = "" if obj is None else str(obj)
            continue
        val = obj.get("value")
        ann = obj.get("annotation")
        if val is not None:
            s = str(val).strip()
            if s or val == 0 or val is False:
                row[str(key)] = s if isinstance(val, str) else str(val)
        if ann is not None and str(ann).strip():
            row[f"{key}__annotation"] = str(ann).strip()

    return row


def collect_fieldnames(rows: Iterable[Dict[str, str]]) -> List[str]:
    keys: Set[str] = set()
    for r in rows:
        keys.update(r.keys())
    rest = sorted(k for k in keys if k not in CORE_ORDER)
    return [k for k in CORE_ORDER if k in keys] + rest


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--output-csv",
        default="/home/mikej10/advbfx/projects/flow_public_samples_pull_v3.csv",
        help="Output CSV path",
    )
    ap.add_argument("--username", default=os.environ.get("FLOWBIO_USERNAME", ""))
    ap.add_argument("--password", default=os.environ.get("FLOWBIO_PASSWORD", ""))
    ap.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    ap.add_argument("--max-projects", type=int, default=0, help="If >0, only first N projects (debug)")
    ap.add_argument(
        "--project-id",
        default="",
        help="If set, pull only this project ID (e.g. 548481478754251364)",
    )
    ap.add_argument(
        "--include-private",
        action="store_true",
        help="Include samples with private=True (default: only public)",
    )
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")

    if not args.username or not args.password:
        logging.error("Set FLOWBIO_USERNAME/FLOWBIO_PASSWORD or use --username/--password")
        return 2

    session = requests.Session()
    token = rest_login(session, args.username, args.password)

    projects = fetch_all_projects(session, token)
    if args.project_id:
        selected = str(args.project_id).strip()
        projects = [p for p in projects if str(p.get("id") or "") == selected]
        if not projects:
            logging.error("project_id %s not found in /projects listing", selected)
            return 2
        logging.info("Using project filter: %s", selected)
    elif args.max_projects > 0:
        projects = projects[: args.max_projects]
    logging.info("Discovered %d project(s) to process", len(projects))

    tasks: List[Tuple[str, str, Optional[bool]]] = []
    for p in projects:
        pid = str(p.get("id") or "")
        if not pid:
            continue
        try:
            slist = fetch_project_samples_list(session, token, pid)
        except requests.HTTPError as e:
            logging.warning("Skip project %s: %s", pid, e)
            continue
        for s in slist:
            sid = s.get("id")
            if not sid:
                continue
            priv = s.get("private")
            tasks.append((pid, str(sid), priv if isinstance(priv, bool) else None))

    seen_sid: Set[str] = set()
    deduped: List[Tuple[str, str, Optional[bool]]] = []
    for pid, sid, p in tasks:
        if sid in seen_sid:
            continue
        seen_sid.add(sid)
        deduped.append((pid, sid, p))
    tasks = deduped
    logging.info("After dedupe: %d unique sample ids", len(tasks))

    if not args.include_private:
        # Drop samples known private from list endpoint when the flag is present
        tasks = [(pid, sid, p) for pid, sid, p in tasks if p is not True]

    rows_out: List[Dict[str, str]] = []
    errors: List[str] = []

    def one_sample(sid: str) -> Tuple[Dict[str, str], Optional[str]]:
        try:
            s2 = requests.Session()
            d = fetch_sample_detail(s2, token, sid)
            if not args.include_private and d.get("private") is True:
                return {}, None
            return flatten_sample_detail(d), None
        except Exception as e:  # noqa: BLE001
            return {}, f"{sid}: {e}"

    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(one_sample, sid): sid for _pid, sid, _p in tasks}
        done = 0
        for fut in as_completed(futs):
            sid = futs[fut]
            row, err = fut.result()
            done += 1
            if done % 500 == 0:
                logging.info("Fetched %d / %d samples", done, len(futs))
            if err:
                errors.append(err)
                if len(errors) <= 25:
                    logging.warning("%s", err)
            elif row:
                rows_out.append(row)

    logging.info("Collected %d public sample rows after detail fetch", len(rows_out))

    if not rows_out:
        logging.error("No rows to write")
        return 3

    fieldnames = collect_fieldnames(rows_out)
    with open(args.output_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in rows_out:
            w.writerow(r)

    manifest = args.output_csv.replace(".csv", "_manifest.json")
    try:
        with open(manifest, "w", encoding="utf-8") as mf:
            json.dump(
                {
                    "row_count": len(rows_out),
                    "project_count": len(projects),
                    "errors": errors[:200],
                    "error_count": len(errors),
                },
                mf,
                indent=2,
            )
    except OSError:
        pass

    logging.info("Wrote %s (%d columns)", args.output_csv, len(fieldnames))
    if errors:
        logging.warning("%d sample fetch errors (first 5 logged at INFO)", len(errors))
        for e in errors[:5]:
            logging.info("%s", e)
    return 0 if not errors else 1


if __name__ == "__main__":
    sys.exit(main())

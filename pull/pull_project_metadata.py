#!/usr/bin/env python3
"""Pull flat metadata CSV for one Flow project (post-upload baseline)."""

from __future__ import annotations

import argparse
import csv
import logging
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from flow_public_samples_pull_v3 import (  # noqa: E402
    collect_fieldnames,
    fetch_sample_detail,
    flatten_sample_detail,
)

DEFAULT_API_BASE = os.environ.get("FLOWBIO_API_BASE", "https://app.flow.bio/api")


def _api_base() -> str:
    return os.environ.get("FLOWBIO_API_BASE", DEFAULT_API_BASE).rstrip("/")


def rest_login_with_base(session: requests.Session, username: str, password: str) -> str:
    r = session.post(
        f"{_api_base()}/login",
        json={"username": username, "password": password},
        timeout=90,
    )
    r.raise_for_status()
    tok = r.json().get("token")
    if not tok:
        raise RuntimeError("Login failed: no token")
    return str(tok)


def fetch_project_samples(
    session: requests.Session, token: str, project_id: str, page_size: int = 100
) -> list[dict]:
    headers = {"Authorization": f"Bearer {token}"}
    page = 1
    out: list[dict] = []
    while True:
        r = session.get(
            f"{_api_base()}/projects/{project_id}/samples",
            params={"page": page, "count": page_size},
            headers=headers,
            timeout=90,
        )
        r.raise_for_status()
        chunk = r.json().get("samples") or []
        if not chunk:
            break
        out.extend(chunk)
        if len(chunk) < page_size:
            break
        page += 1
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project-id", required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--username", default=os.environ.get("FLOWBIO_USERNAME", ""))
    ap.add_argument("--password", default=os.environ.get("FLOWBIO_PASSWORD", ""))
    ap.add_argument("--include-private", action="store_true")
    ap.add_argument("--workers", type=int, default=12)
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")

    if not args.username or not args.password:
        logging.error("Set FLOWBIO_USERNAME/FLOWBIO_PASSWORD")
        return 2

    session = requests.Session()
    token = rest_login_with_base(session, args.username, args.password)
    samples = fetch_project_samples(session, token, args.project_id)
    if not args.include_private:
        samples = [s for s in samples if s.get("private") is not True]
    logging.info("Project %s: %d sample(s) to fetch", args.project_id, len(samples))

    rows: list[dict[str, str]] = []
    errors: list[str] = []

    def one(sid: str) -> tuple[dict[str, str], str | None]:
        try:
            s2 = requests.Session()
            import flow_public_samples_pull_v3 as pull_mod

            old = pull_mod.API_BASE
            pull_mod.API_BASE = _api_base()
            try:
                d = fetch_sample_detail(s2, token, sid)
            finally:
                pull_mod.API_BASE = old
            if not args.include_private and d.get("private") is True:
                return {}, None
            return flatten_sample_detail(d), None
        except Exception as exc:  # noqa: BLE001
            return {}, f"{sid}: {exc}"

    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(one, str(s["id"])): s for s in samples if s.get("id")}
        for fut in as_completed(futs):
            row, err = fut.result()
            if err:
                errors.append(err)
            elif row:
                rows.append(row)

    rows.sort(key=lambda r: r.get("sample_name", ""))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = collect_fieldnames(rows)
    with args.output.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    logging.info("Wrote %s (%d rows, %d errors)", args.output, len(rows), len(errors))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())

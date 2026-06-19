#!/usr/bin/env python3
import argparse
import getpass
import json
import logging
import re
from typing import Any, Dict, List, Optional, Pattern, Set

import requests

API_BASE = "https://api.flow.bio"
APP_API_BASE = "https://app.flow.bio/api"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Log in to Flow.bio, fetch all samples in a project, load full sample details, "
            "and collect sample-associated data metadata."
        )
    )
    parser.add_argument(
        "-p",
        "--project",
        required=True,
        help="Flow project URL or project ID (e.g. https://app.flow.bio/projects/<id>/ or <id>)",
    )
    parser.add_argument("-u", "--username", default=None, help="Flow username (prompted if omitted)")
    parser.add_argument("--password", default=None, help="Flow password (prompted if omitted)")
    parser.add_argument("--page-size", type=int, default=100, help="Page size for project sample listing")
    parser.add_argument(
        "--data-page-size",
        type=int,
        default=100,
        help="Page size for /samples/{id}/data listing (default: 100)",
    )
    parser.add_argument(
        "--max-samples",
        type=int,
        default=None,
        help="Optional cap on number of project samples to process",
    )
    parser.add_argument(
        "--uploaded-only",
        action="store_true",
        help="Only emit samples that have at least one associated data record marked as Uploaded",
    )
    parser.add_argument(
        "--data-source",
        choices=("filesets", "sample_data", "both"),
        default="filesets",
        help=(
            "Where to list sample-associated files: filesets (uploaded FASTQs), "
            "sample_data (/samples/{id}/data, includes pipeline outputs), or both."
        ),
    )
    parser.add_argument(
        "--filename-regex",
        default=None,
        help=(
            "Only keep data files whose filename matches this regex. "
            "Example for C34_2.markdup.sorted.bam: '\\.markdup\\.sorted\\.bam$' "
            "(in shell: --filename-regex '\\.markdup\\.sorted\\.bam$')"
        ),
    )
    parser.add_argument(
        "--output-json",
        default=None,
        help="Optional path to save collected sample/data metadata as JSON",
    )
    parser.add_argument("--debug", action="store_true", help="Enable debug logging")
    return parser.parse_args()


def setup_logging(debug: bool = False) -> None:
    level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s | %(levelname)-8s | %(message)s")


def extract_project_id(project_input: str) -> str:
    candidate = project_input.strip()
    if candidate.isdigit():
        return candidate
    match = re.search(r"/projects/(\d+)", candidate)
    if match:
        return match.group(1)
    raise ValueError(f"Could not parse project id from input: {project_input}")


def raise_for_status(resp: requests.Response) -> None:
    try:
        resp.raise_for_status()
    except requests.HTTPError as exc:
        message = f"HTTP {resp.status_code} error for {resp.request.method} {resp.url}: {resp.text}"
        raise requests.HTTPError(message) from exc


def rest_login(session: requests.Session, username: Optional[str], password: Optional[str]) -> str:
    if not username:
        username = input("Enter your username: ").strip()
    if not password:
        password = getpass.getpass("Enter your password: ")
    resp = session.post(
        f"{API_BASE}/login",
        json={"username": username, "password": password},
        timeout=30,
    )
    raise_for_status(resp)
    payload = resp.json()
    token = payload.get("token")
    if not token:
        raise RuntimeError("Login succeeded but no token returned in response body")
    return token


def fetch_all_project_samples(
    session: requests.Session, token: str, project_id: str, page_size: int = 100
) -> List[Dict[str, Any]]:
    headers = {"Authorization": f"Bearer {token}"}
    page = 1
    collected: List[Dict[str, Any]] = []
    while True:
        resp = session.get(
            f"{API_BASE}/projects/{project_id}/samples",
            params={"page": page, "count": page_size},
            headers=headers,
            timeout=30,
        )
        raise_for_status(resp)
        payload = resp.json()
        samples = payload.get("samples", [])
        if not samples:
            break
        collected.extend(samples)
        logging.info("Fetched project sample page %s (%s samples)", page, len(samples))
        if len(samples) < page_size:
            break
        page += 1
    logging.info("Collected %s project sample records", len(collected))
    return collected


def fetch_sample_detail(session: requests.Session, token: str, sample_id: str) -> Dict[str, Any]:
    headers = {"Authorization": f"Bearer {token}"}
    resp = session.get(f"{API_BASE}/samples/{sample_id}", headers=headers, timeout=30)
    raise_for_status(resp)
    return resp.json()


def fetch_data_detail(session: requests.Session, token: str, data_id: str) -> Dict[str, Any]:
    headers = {"Authorization": f"Bearer {token}"}
    resp = session.get(f"{API_BASE}/data/{data_id}", headers=headers, timeout=30)
    raise_for_status(resp)
    return resp.json()


def _collect_sample_data_entries(
    session: requests.Session, token: str, sample_id: str, page_size: int = 100
) -> List[Dict[str, Any]]:
    headers = {"Authorization": f"Bearer {token}"}
    page = 1
    collected: List[Dict[str, Any]] = []
    seen_ids: Set[str] = set()
    total_reported: Optional[int] = None

    while True:
        resp = session.get(
            f"{APP_API_BASE}/samples/{sample_id}/data",
            params={"page": page, "count": page_size},
            headers=headers,
            timeout=30,
        )
        raise_for_status(resp)
        payload = resp.json()
        if total_reported is None:
            total_reported = payload.get("count")
        entries = payload.get("data", [])
        if not isinstance(entries, list) or not entries:
            break

        for entry in entries:
            if not isinstance(entry, dict):
                continue
            data_id = entry.get("id")
            if not data_id:
                continue
            sid = str(data_id)
            if sid in seen_ids:
                continue
            seen_ids.add(sid)
            collected.append(entry)

        if len(entries) < page_size:
            break
        if total_reported is not None and len(collected) >= int(total_reported):
            break
        page += 1

    logging.debug(
        "Sample %s: paginated data entries=%s (reported total=%s, pages=%s)",
        sample_id,
        len(collected),
        total_reported,
        page,
    )
    return collected


def _collect_data_ids_from_filesets(sample_detail: Dict[str, Any]) -> List[str]:
    data_ids: List[str] = []
    filesets = sample_detail.get("filesets", [])
    if not isinstance(filesets, list):
        return data_ids
    for fileset in filesets:
        if not isinstance(fileset, dict):
            continue
        data_entries = fileset.get("data", [])
        if not isinstance(data_entries, list):
            continue
        for entry in data_entries:
            if not isinstance(entry, dict):
                continue
            data_id = entry.get("id")
            if data_id:
                data_ids.append(str(data_id))
    return data_ids


def _is_uploaded_data_record(record: Dict[str, Any]) -> bool:
    if record.get("execution") is None:
        return True
    possible_keys = ("source", "origin", "status", "data_origin", "upload_source", "type")
    for key in possible_keys:
        value = record.get(key)
        if isinstance(value, str) and value.strip().lower() == "uploaded":
            return True
    provenance = record.get("provenance")
    if isinstance(provenance, dict):
        for key in ("source", "origin", "status"):
            value = provenance.get(key)
            if isinstance(value, str) and value.strip().lower() == "uploaded":
                return True
    return False


def _matches_filename(record: Dict[str, Any], pattern: Optional[Pattern[str]]) -> bool:
    if pattern is None:
        return True
    filename = record.get("filename") or ""
    return bool(pattern.search(str(filename)))


def collect_sample_and_data_metadata(
    project_samples: List[Dict[str, Any]],
    session: requests.Session,
    token: str,
    uploaded_only: bool = False,
    data_source: str = "filesets",
    filename_pattern: Optional[Pattern[str]] = None,
    data_page_size: int = 100,
    max_samples: Optional[int] = None,
) -> List[Dict[str, Any]]:
    collected: List[Dict[str, Any]] = []
    total = len(project_samples) if max_samples is None else min(len(project_samples), max_samples)
    for idx, sample in enumerate(project_samples, start=1):
        if max_samples is not None and idx > max_samples:
            break
        sample_id = str(sample.get("id", ""))
        sample_name = sample.get("name", "")
        if not sample_id:
            continue
        detail = fetch_sample_detail(session, token, sample_id)
        sample_metadata = detail.get("metadata") if isinstance(detail.get("metadata"), dict) else {}
        data_ids: List[str] = []
        if data_source in ("filesets", "both"):
            data_ids.extend(_collect_data_ids_from_filesets(detail))
        list_entries: List[Dict[str, Any]] = []
        if data_source in ("sample_data", "both"):
            try:
                list_entries = _collect_sample_data_entries(
                    session, token, sample_id, page_size=data_page_size
                )
                data_ids.extend(str(entry["id"]) for entry in list_entries if entry.get("id"))
            except Exception as exc:
                logging.warning("Failed to fetch /samples/%s/data: %s", sample_id, exc)
        unique_ids: List[str] = []
        seen = set()
        for data_id in data_ids:
            if data_id not in seen:
                seen.add(data_id)
                unique_ids.append(data_id)

        ids_to_fetch: List[str] = []
        if filename_pattern is not None and list_entries:
            by_filename: Dict[str, Dict[str, Any]] = {}
            for entry in list_entries:
                if not _matches_filename(entry, filename_pattern):
                    continue
                filename = str(entry.get("filename") or "")
                if not filename:
                    continue
                prev = by_filename.get(filename)
                if prev is None or int(entry.get("size") or 0) >= int(prev.get("size") or 0):
                    by_filename[filename] = entry
            ids_to_fetch = [str(entry["id"]) for entry in by_filename.values() if entry.get("id")]
        else:
            ids_to_fetch = unique_ids

        data_candidates: List[Dict[str, Any]] = []
        for data_id in ids_to_fetch:
            try:
                data_candidates.append(fetch_data_detail(session, token, data_id))
            except Exception as exc:
                logging.warning("Failed to fetch /data/%s: %s", data_id, exc)
        if filename_pattern is not None and not list_entries:
            data_candidates = [d for d in data_candidates if _matches_filename(d, filename_pattern)]
        uploaded_data = [d for d in data_candidates if _is_uploaded_data_record(d)]

        if uploaded_only and not uploaded_data:
            logging.debug("Skipping sample %s with no Uploaded data records", sample_id)
            continue
        if filename_pattern is not None and not data_candidates:
            logging.debug("Skipping sample %s with no filename matches", sample_id)
            continue

        record = {
            "sample_id": sample_id,
            "sample_name": sample_name or detail.get("name", ""),
            "sample_metadata": sample_metadata,
            "data_count": len(data_candidates),
            "uploaded_data_count": len(uploaded_data),
            "data": data_candidates,
            "uploaded_data": uploaded_data,
        }
        collected.append(record)
        logging.info(
            "Processed sample %s/%s: id=%s, data=%s, uploaded=%s",
            idx,
            total,
            sample_id,
            len(data_candidates),
            len(uploaded_data),
        )
    return collected


def main() -> None:
    args = parse_args()
    setup_logging(args.debug)
    project_id = extract_project_id(args.project)

    filename_pattern = re.compile(args.filename_regex) if args.filename_regex else None

    session = requests.Session()
    token = rest_login(session, args.username, args.password)
    project_samples = fetch_all_project_samples(
        session=session, token=token, project_id=project_id, page_size=args.page_size
    )
    collected = collect_sample_and_data_metadata(
        project_samples=project_samples,
        session=session,
        token=token,
        uploaded_only=args.uploaded_only,
        data_source=args.data_source,
        filename_pattern=filename_pattern,
        data_page_size=args.data_page_size,
        max_samples=args.max_samples,
    )

    print(f"\nProject {project_id}: collected {len(collected)} sample detail records")
    if collected:
        first = collected[0]
        print(
            "First sample preview: "
            f"id={first['sample_id']} name={first['sample_name']} "
            f"data_count={first['data_count']} uploaded_data_count={first['uploaded_data_count']}"
        )
        if first["uploaded_data_count"] > 0:
            uploaded_ids = [x.get("id") for x in first["uploaded_data"][:3]]
            print(f"First sample uploaded data IDs (up to 3): {uploaded_ids}")

    if args.output_json:
        with open(args.output_json, "w", encoding="utf-8") as handle:
            json.dump(collected, handle, indent=2)
        logging.info("Wrote JSON output to %s", args.output_json)


if __name__ == "__main__":
    main()

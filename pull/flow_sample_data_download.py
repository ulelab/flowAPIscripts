#!/usr/bin/env python3
import argparse
import getpass
import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

import requests

API_BASE_DEFAULT = "https://api.flow.bio"
APP_API_BASE_DEFAULT = "https://app.flow.bio/api"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Bulk download Flow data IDs from a text file. "
            "Reads one data ID per line, logs in, fetches metadata, and downloads files."
        )
    )
    parser.add_argument(
        "-i",
        "--ids-file",
        required=True,
        help="Path to data_id.txt (one data_id per line; # comments allowed)",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        default="flow_downloads",
        help="Output directory for downloaded files",
    )
    parser.add_argument("-u", "--username", default=None, help="Flow username (prompted if omitted)")
    parser.add_argument("--password", default=None, help="Flow password (prompted if omitted)")
    parser.add_argument(
        "--api-base",
        default=API_BASE_DEFAULT,
        help=f"API base URL (default: {API_BASE_DEFAULT})",
    )
    parser.add_argument(
        "--app-api-base",
        default=APP_API_BASE_DEFAULT,
        help=f"App API base URL (default: {APP_API_BASE_DEFAULT})",
    )
    parser.add_argument(
        "--download-endpoint-template",
        default=None,
        help=(
            "Override download URL template, e.g. 'https://app.flow.bio/api/data/{data_id}/download'. "
            "Use {data_id} placeholder."
        ),
    )
    parser.add_argument(
        "--manifest-json",
        default=None,
        help="Optional path to write metadata manifest JSON",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Do not download files; only validate IDs and print/write manifest",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing files in outdir",
    )
    parser.add_argument("--debug", action="store_true", help="Enable debug logging")
    return parser.parse_args()


def setup_logging(debug: bool) -> None:
    level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s | %(levelname)-8s | %(message)s")


def read_data_ids(ids_file: Path) -> List[str]:
    if not ids_file.exists():
        raise FileNotFoundError(f"IDs file not found: {ids_file}")
    seen = set()
    data_ids: List[str] = []
    for raw_line in ids_file.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        # Accept "data_id<TAB>anything", keep first token.
        data_id = line.split()[0]
        if data_id not in seen:
            seen.add(data_id)
            data_ids.append(data_id)
    if not data_ids:
        raise ValueError(f"No data IDs found in file: {ids_file}")
    return data_ids


def raise_for_status(resp: requests.Response) -> None:
    try:
        resp.raise_for_status()
    except requests.HTTPError as exc:
        msg = f"HTTP {resp.status_code} error for {resp.request.method} {resp.url}: {resp.text[:500]}"
        raise requests.HTTPError(msg) from exc


def rest_login(session: requests.Session, api_base: str, username: Optional[str], password: Optional[str]) -> str:
    if not username:
        username = input("Enter your username: ").strip()
    if not password:
        password = getpass.getpass("Enter your password: ")
    resp = session.post(
        f"{api_base}/login",
        json={"username": username, "password": password},
        timeout=30,
    )
    raise_for_status(resp)
    payload = resp.json()
    token = payload.get("token")
    if not token:
        raise RuntimeError("Login response did not contain token")
    return token


def fetch_data_metadata(session: requests.Session, api_base: str, token: str, data_id: str) -> Dict[str, Any]:
    headers = {"Authorization": f"Bearer {token}"}
    resp = session.get(f"{api_base}/data/{data_id}", headers=headers, timeout=30)
    raise_for_status(resp)
    return resp.json()


def _default_download_urls(api_base: str, app_api_base: str, data_id: str) -> List[str]:
    return [
        f"https://app.flow.bio/files/downloads/{data_id}",
        f"{api_base}/data/{data_id}/download",
        f"{app_api_base}/data/{data_id}/download",
        f"{api_base}/data/{data_id}/content",
    ]


def _looks_like_file_response(resp: requests.Response) -> bool:
    ctype = (resp.headers.get("Content-Type") or "").lower()
    if "application/json" in ctype or "text/html" in ctype:
        return False
    return True


def stream_download(
    session: requests.Session,
    token: str,
    metadata: Dict[str, Any],
    outdir: Path,
    api_base: str,
    app_api_base: str,
    endpoint_template: Optional[str],
    overwrite: bool,
) -> Dict[str, Any]:
    data_id = str(metadata["id"])
    filename = metadata.get("filename") or f"{data_id}.bin"
    destination = outdir / filename
    if destination.exists() and not overwrite:
        return {"status": "skipped_exists", "data_id": data_id, "path": str(destination)}

    headers = {"Authorization": f"Bearer {token}"}
    if endpoint_template:
        urls = [endpoint_template.format(data_id=data_id)]
    else:
        urls = _default_download_urls(api_base=api_base, app_api_base=app_api_base, data_id=data_id)
        # Flow file downloads often require filename in the URL path.
        if metadata.get("filename"):
            urls.insert(0, f"https://app.flow.bio/files/downloads/{data_id}/{metadata['filename']}")

    last_error = None
    for url in urls:
        try:
            with session.get(url, headers=headers, stream=True, timeout=120, allow_redirects=True) as resp:
                if resp.status_code != 200:
                    last_error = f"HTTP {resp.status_code} from {url}"
                    continue
                if not _looks_like_file_response(resp):
                    snippet = resp.text[:200]
                    last_error = f"Non-file response from {url}: content-type={resp.headers.get('Content-Type')} body={snippet}"
                    continue
                outdir.mkdir(parents=True, exist_ok=True)
                with open(destination, "wb") as handle:
                    for chunk in resp.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            handle.write(chunk)
                return {"status": "downloaded", "data_id": data_id, "path": str(destination), "url": url}
        except Exception as exc:
            last_error = f"{type(exc).__name__}: {exc}"
            continue

    return {
        "status": "failed",
        "data_id": data_id,
        "filename": filename,
        "error": last_error or "No usable download endpoint produced a file response",
        "attempted_urls": urls,
    }


def main() -> None:
    args = parse_args()
    setup_logging(args.debug)

    ids_file = Path(args.ids_file)
    outdir = Path(args.outdir)
    data_ids = read_data_ids(ids_file)
    logging.info("Loaded %s unique data IDs from %s", len(data_ids), ids_file)

    session = requests.Session()
    token = rest_login(session, args.api_base, args.username, args.password)

    manifest: List[Dict[str, Any]] = []
    for idx, data_id in enumerate(data_ids, start=1):
        try:
            md = fetch_data_metadata(session, args.api_base, token, data_id)
            record = {
                "data_id": data_id,
                "filename": md.get("filename"),
                "size": md.get("size"),
                "sample_id": (md.get("sample") or {}).get("id") if isinstance(md.get("sample"), dict) else None,
                "sample_name": (md.get("sample") or {}).get("name") if isinstance(md.get("sample"), dict) else None,
                "project_id": (md.get("project") or {}).get("id") if isinstance(md.get("project"), dict) else None,
                "project_name": (md.get("project") or {}).get("name") if isinstance(md.get("project"), dict) else None,
                "execution": md.get("execution"),
                "absolute_path": md.get("absolute_path"),
            }
            if args.dry_run:
                record["download_status"] = "dry_run"
            else:
                result = stream_download(
                    session=session,
                    token=token,
                    metadata=md,
                    outdir=outdir,
                    api_base=args.api_base,
                    app_api_base=args.app_api_base,
                    endpoint_template=args.download_endpoint_template,
                    overwrite=args.overwrite,
                )
                record["download_status"] = result.get("status")
                record["download_result"] = result
            manifest.append(record)
            logging.info("[%s/%s] data_id=%s filename=%s status=%s", idx, len(data_ids), data_id, record["filename"], record["download_status"])
        except Exception as exc:
            logging.error("[%s/%s] data_id=%s failed metadata/download: %s", idx, len(data_ids), data_id, exc)
            manifest.append({"data_id": data_id, "download_status": "failed", "error": str(exc)})

    if args.manifest_json:
        manifest_path = Path(args.manifest_json)
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        logging.info("Wrote manifest to %s", manifest_path)

    ok = sum(1 for x in manifest if x.get("download_status") in {"downloaded", "dry_run", "skipped_exists"})
    fail = sum(1 for x in manifest if x.get("download_status") == "failed")
    print(f"\nCompleted {len(manifest)} IDs: ok={ok}, failed={fail}")
    if fail:
        print("Failures were recorded in manifest output (or log) with endpoint attempts.")


if __name__ == "__main__":
    main()

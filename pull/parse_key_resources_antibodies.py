#!/usr/bin/env python3
"""Parse Mol Cell Key Resources antibody rows into purification_agent proposals TSV."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

# GSE290281 eCLIP IP targets (INPUT/SMInput omitted — purification_agent stays empty)
DEFAULT_TARGETS = [
    "CPSF5",
    "CPSF6",
    "EIF4B",
    "GRB2",
    "LGALS3",
    "MBNL2",
    "RBM10",
    "RBM22",
    "RNPS1",
]

CATALOG_RE = re.compile(r"\b(A\d{3}-\d{3}[A-Z]|[A-Z]{1,2}\d{3,6}[A-Z]?)\b")
TARGET_ALIASES = {
    "NUDT21": "CPSF5",
    "CFIM25": "CPSF5",
    "CFIM68": "CPSF6",
}


def _norm_species(token: str) -> str:
    t = token.strip().lower()
    if t.startswith("rabbit"):
        return "Rabbit"
    if t.startswith("mouse"):
        return "Mouse"
    if t.startswith("goat"):
        return "Goat"
    return token.strip().title()


def format_agent(species: str, target: str, vendor: str, catalog: str) -> str:
    vendor = vendor.strip()
    catalog = catalog.strip()
    cat_part = f"{vendor} {catalog}".strip() if vendor else catalog
    return f"{species} Anti-{target} ({cat_part})"


def parse_line(line: str) -> tuple[str, str, str] | None:
    """Return (target, agent_string, evidence_line) or None."""
    raw = line.strip()
    if not raw or raw.lower().startswith("reagent"):
        return None
    lower = raw.lower()
    if "antibod" not in lower and "anti-" not in lower and "anti " not in lower:
        return None

    catalogs = CATALOG_RE.findall(raw)
    if not catalogs:
        return None
    catalog = catalogs[0]

    target = ""
    for alias, canonical in TARGET_ALIASES.items():
        if re.search(rf"\b{re.escape(alias)}\b", raw, re.I):
            target = canonical
            break
    if not target:
        for t in DEFAULT_TARGETS:
            if re.search(rf"\b{re.escape(t)}\b", raw, re.I):
                target = t
                break
    if not target:
        return None

    species = "Rabbit"
    for sp in ("rabbit", "mouse", "goat"):
        if sp in lower:
            species = _norm_species(sp)
            break

    vendor = ""
    for v in (
        "Bethyl Laboratories",
        "Bethyl",
        "Sigma-Aldrich",
        "Sigma",
        "Abcam",
        "Proteintech",
        "Cell Signaling",
        "Santa Cruz",
        "Thermo Fisher",
    ):
        if v.lower() in lower:
            vendor = "Sigma" if v.startswith("Sigma") else v.replace(" Laboratories", "")
            break

    agent = format_agent(species, target, vendor, catalog)
    return target, agent, raw


def parse_text(text: str) -> dict[str, tuple[str, str]]:
    """target -> (agent, evidence). Last match wins."""
    out: dict[str, tuple[str, str]] = {}
    for line in text.splitlines():
        parsed = parse_line(line)
        if parsed:
            target, agent, evidence = parsed
            out[target] = (agent, evidence)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "input",
        type=Path,
        help="Pasted Key resources / antibody table (plain text)",
    )
    ap.add_argument(
        "--output",
        type=Path,
        required=True,
        help="purification_agent_proposals.tsv",
    )
    ap.add_argument(
        "--status",
        default="pending_confirmation",
        choices=("pending_confirmation", "confirmed"),
    )
    args = ap.parse_args()

    found = parse_text(args.input.read_text(encoding="utf-8", errors="replace"))
    rows = []
    for target in DEFAULT_TARGETS:
        if target in found:
            agent, evidence = found[target]
            rows.append(
                {
                    "purification_target": target,
                    "proposed_purification_agent": agent,
                    "evidence_quote": evidence[:500],
                    "status": args.status,
                }
            )
        else:
            rows.append(
                {
                    "purification_target": target,
                    "proposed_purification_agent": "",
                    "evidence_quote": "NOT FOUND in input — check Key resources table",
                    "status": "pending_confirmation",
                }
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "purification_target",
                "proposed_purification_agent",
                "evidence_quote",
                "status",
            ],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(rows)

    missing = [r["purification_target"] for r in rows if not r["proposed_purification_agent"]]
    print(f"Wrote {args.output} ({len(rows) - len(missing)}/{len(rows)} targets matched)")
    if missing:
        print("Missing:", ", ".join(missing))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

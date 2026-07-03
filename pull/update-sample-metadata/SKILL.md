---
name: update-sample-metadata
description: >-
  Pull live metadata for uploaded Flow.bio samples, draft purification-agent
  (and other field) updates from paper Methods, and push changes via the v2
  metadata push script. Use when samples are already on Flow and need corrected
  Purification Agent, purification target annotation, or other whitelisted fields
  after upload; when GEO only says "RBP-specific antibody"; or when the user
  asks to patch sample metadata post-upload.
---

# Update Sample Metadata (post-upload)

## When to use

- Samples are **already uploaded**; annotation at upload time was incomplete (common for **Purification Agent**).
- The specific antibody (vendor + catalog) is in the **paper Methods** (often a table under eCLIP / CLIP), not in GEO.
- Same pattern as `purification_target__annotation`: agent reads source text → human confirms → script pushes.

## Purification Agent vs Purification Target Annotation

| Field | API key | Transport | Meaning |
|-------|---------|-----------|---------|
| Purification Target | `purification_target` | GraphQL / REST parent | Gene symbol (`CPSF5`, `GRB2`) or `SMInput` for size-matched input |
| Purification Target Annotation | `purification_target__annotation` | REST `/edit` | Tag on construct (`cV5`, `3xFLAG`) |
| **Purification Agent** | `purification_agent` | **GraphQL** | Full antibody descriptor |

**Agent job:** extract antibody strings from Methods (CLIP subsection). Do **not** invent catalog numbers.

Formatting rules: [reference/purification-agent-format.md](reference/purification-agent-format.md)

## Workflow

```
1. Pull baseline CSV     →  samples_baseline.csv
2. Agent reads Methods   →  purification_agent_proposals.tsv (or field-specific proposals)
3. Human confirms        →  status=confirmed on each row
4. Build updated CSV     →  samples_updated.csv
5. Dry-run diff          →  flow_public_samples_push_metadata_v2.py --dry-run
6. Push + verify         →  --yes (or interactive)
```

### Step 1 — Pull project samples

```bash
export FLOWBIO_USERNAME=... FLOWBIO_PASSWORD=...

python3 flowAPIscripts/pull/pull_project_metadata.py \
  --project-id 838933490352991054 \
  --output projects/GSE290281/flow-output/samples_baseline.csv \
  --include-private
```

Uses `flatten_sample_detail` from `flow_public_samples_pull_v3.py`. Defaults to `https://app.flow.bio/api` (set `FLOWBIO_API_BASE` if needed).

Optional audit of uploaded FASTQs:

```bash
python3 flowAPIscripts/pull/flow_sample_data_pull.py \
  -p 838933490352991054 \
  --uploaded-only \
  --output-json projects/GSE290281/flow-output/samples_data_audit.json
```

### Step 2 — Agent: Methods / Key resources table → proposals

1. Resolve the linked publication from `flagged_papers.json`, GEO `!Series_pubmed_id`, or project `pubmed` column.
2. Open **STAR Methods → Key resources table** (Mol Cell). For PMID **42361791** / [ScienceDirect](https://www.sciencedirect.com/science/article/pii/S1097276526003825) — antibody catalog numbers are in that table, not PubMed abstract.
3. Copy antibody rows (species, target, vendor, catalog) into `paper_key_resources_antibodies.txt`.
4. Parse into proposals:

```bash
python3 flowAPIscripts/pull/parse_key_resources_antibodies.py \
  projects/GSE290281/flow-output/paper_key_resources_antibodies.txt \
  --output projects/GSE290281/flow-output/purification_agent_proposals.tsv
```

5. Review TSV; set `status` to `confirmed` on each IP target row.
6. **INPUT / SMInput** samples: leave `purification_agent` empty (no row needed — `apply_metadata_proposals.py` skips `SMInput`).

| column | required | notes |
|--------|----------|-------|
| `purification_target` | yes | e.g. `CPSF5`, `GRB2` |
| `proposed_purification_agent` | for IP | formatted string; empty for SMInput |
| `evidence_quote` | yes | short Methods quote |
| `status` | yes | `pending_confirmation` → `confirmed` |

Per-sample overrides: add `sample_id` column (see script `--by-sample`).

Save paper excerpt optionally as `paper_methods_antibodies.txt` for audit.

### Step 3 — Build updated CSV

```bash
python3 flowAPIscripts/pull/apply_metadata_proposals.py \
  --baseline projects/GSE290281/flow-output/samples_baseline.csv \
  --proposals projects/GSE290281/flow-output/purification_agent_proposals.tsv \
  --output projects/GSE290281/flow-output/samples_updated.csv \
  --field purification_agent
```

Writes `CONFIRM_METADATA_UPDATES.md` listing rows that will change.

### Step 4 — Dry-run and push

```bash
python3 flowAPIscripts/pull/flow_public_samples_push_metadata_v2.py \
  --dry-run \
  --baseline projects/GSE290281/flow-output/samples_baseline.csv \
  --updated projects/GSE290281/flow-output/samples_updated.csv

python3 flowAPIscripts/pull/flow_public_samples_push_metadata_v2.py \
  --yes \
  --baseline projects/GSE290281/flow-output/samples_baseline.csv \
  --updated projects/GSE290281/flow-output/samples_updated.csv \
  --audit-dir projects/GSE290281/flow-output/metadata_push_audit
```

`purification_agent` updates use **GraphQL** `updateSample` (`purificationAgent`). Annotation fields use REST `POST /samples/{id}/edit`.

### Step 5 — Verify

- Push script re-fetches each sample and compares live vs expected.
- Re-pull baseline and spot-check IP rows in Flow UI.

## Other pushable fields

From `flow_public_samples_push_metadata_v2.py`:

- **GraphQL:** `sample_name`, `condition`, `comments`, `experimental_method`, `purification_agent`, `purification_target`, `source`
- **REST `/edit`:** `purification_target__annotation`, `source__annotation`

Use the same baseline → proposals → updated → push pattern; extend `apply_metadata_proposals.py` `--field` or add columns to proposals TSV.

## Guardrails

- Never push without `--dry-run` review when more than one field changes.
- Do not clear existing metadata unless `--allow-clear` is intentional.
- `sample_id` in baseline CSV is the join key — do not rename samples in the same pass unless intended.
- Credentials: `FLOWBIO_USERNAME` / `FLOWBIO_PASSWORD` (never commit).

## References

- [reference/api-and-scripts.md](reference/api-and-scripts.md) — script map and API notes
- [reference/purification-agent-format.md](reference/purification-agent-format.md) — antibody string format
- `flow-compile` barcode agent search (Methods CLIP focus) — same agentic pattern for pre-upload gaps

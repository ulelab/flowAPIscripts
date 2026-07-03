# Flow API scripts for metadata pull / push

## Pull — sample IDs and current metadata

### `flow_public_samples_pull_v3.py` (primary)

- **Path:** `flowAPIscripts/pull/flow_public_samples_pull_v3.py`
- **Login:** `POST https://api.flow.bio/login` (may 503; skill `pull_project_metadata.py` defaults to `https://app.flow.bio/api`)
- **List:** `GET /projects/{id}/samples?page=&count=`
- **Detail:** `GET /samples/{sample_id}`
- **Flatten:** `flatten_sample_detail()` → CSV row with `sample_id`, `purification_agent`, `purification_target`, `purification_target__annotation`, etc.

```bash
python3 flowAPIscripts/pull/flow_public_samples_pull_v3.py \
  --project-id PROJECT_ID \
  --include-private \
  --output-csv samples_baseline.csv
```

### `flow_sample_data_pull.py` (uploaded files audit)

- **Path:** `flowAPIscripts/pull/flow_sample_data_pull.py`
- **Use:** confirm uploads (`--uploaded-only`), list FASTQ filenames per sample, optional `--filename-regex`
- **Does not push metadata** — complements baseline pull when verifying the right files are on Flow.

```bash
python3 flowAPIscripts/pull/flow_sample_data_pull.py \
  -p PROJECT_ID \
  --uploaded-only \
  --output-json samples_data_audit.json
```

### Skill helper

`scripts/pull_project_metadata.py` — thin wrapper: project-scoped pull, `app.flow.bio/api` default, reuses `flatten_sample_detail`.

## Push — diff baseline vs updated CSV

### `flow_public_samples_push_metadata_v2.py`

- **Path:** `flowAPIscripts/pull/flow_public_samples_push_metadata_v2.py`
- **Join key:** `sample_id`
- **Diff:** `compute_pending_changes(baseline, updated)` → per-field `FieldChange` with transport `graphql` or `rest`
- **`purification_agent`:** GraphQL `updateSample(purificationAgent: ...)`
- **`purification_target__annotation`:** REST `POST https://app.flow.bio/api/samples/{id}/edit` with parent + annotation keys

```bash
# Always dry-run first
python3 flow_public_samples_push_metadata_v2.py \
  --dry-run \
  --baseline samples_baseline.csv \
  --updated samples_updated.csv

python3 flow_public_samples_push_metadata_v2.py \
  --yes \
  --baseline samples_baseline.csv \
  --updated samples_updated.csv \
  --audit-dir ./metadata_push_audit
```

Single-field REST test:

```bash
python3 flow_public_samples_push_metadata_v2.py \
  --test-sample-id SAMPLE_ID \
  --test-field purification_target__annotation \
  --test-annotation cV5
```

## CSV requirements

- Both files must contain the same `sample_id` keys for rows you intend to update.
- Updated rows not in baseline are **skipped** with a warning.
- Only whitelisted columns generate pushes; other columns are ignored.

## Credentials

- Env: `FLOWBIO_USERNAME`, `FLOWBIO_PASSWORD`
- Or `--username` / `--password`
- Project `.flow_credentials.env` pattern: `source projects/.../.flow_credentials.env`

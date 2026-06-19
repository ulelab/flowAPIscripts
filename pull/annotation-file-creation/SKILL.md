---
name: annotation-file-creation
description: Convert metadata from CLIP and RNA-Seq studies into the upload annotation table used by uploadsample_flowbio_v6.py and uploadsample_multiple_v6.sbatch. Use when the user asks to create, normalize, or validate annotation spreadsheets, SraRunTable mappings, sample naming, CLIP/RNA-Seq field mapping, GEO/PubMed metadata extraction, strandedness assignment, or paired-end file column handling.
---

# Annotation File Creation

## Purpose
Create and validate annotation tables for FlowBio uploads by mapping study metadata to the template columns used by:
- `flowAPIscripts/upload/uploadsample_flowbio_v6.py`
- `flowAPIscripts/upload/uploadsample_multiple_v6.sbatch`

Use the column conventions from:
- `flowAPIscripts/test-datasets/Testtemplate.xlsx`

## Workflow
1. Determine whether samples are RNA-Seq or CLIP.
2. Identify single-end vs paired-end and map file names accordingly.
3. Build `Sample Name` using the naming rules below.
4. Fill required study metadata fields (project, scientist, PI, organization, method/type, condition, sequencer, GEO/PubMed).
5. Fill species/source fields with normalized values.
6. For CLIP, capture barcode and protein target fields.
7. For RNA-Seq, set strandedness (`forward`, `reverse`, `auto`, or `unstranded`).
8. Leave nonessential columns empty unless the user explicitly requests them.

## Column Mapping Rules

### File columns
- **File name**: FASTQ filename being processed.
- **File 2**: only for paired-end datasets; use mate pair filename.

### Sample Name
- Build sample names in this order:
  1. `Protein target` (CLIP only, if known)
  2. `Cell or tissue type`
  3. `Species`
  4. `Condition(s)`
  5. `Replicate` / other distinguishing identifiers
- Keep names concise but uniquely identifying.

### Study identity fields
- **Project Name**: study name.
- **Scientist**: uploader.
- **PI**: lab head (often last author on paper).
- **Organisation**: institution where the study was conducted.

### Assay description fields
- **Purification agent**: antibody or bead used for IP/purification.
- **Experimental method**: specific protocol (for example `Quant-Seq`, `iCLIP`).
- **Type**: broad modality (`RNA-Seq`, `CLIP`).
- **Condition**: sample-specific controlled variable(s), not the entire protocol narrative.

### Sequencing/platform fields
- **Sequencer**: sequencing machine/platform, often Illumina model.
- Treat values labeled as `Platform` in source metadata as sequencer input.

### Barcode / primer fields
- Most barcode/primer fields are optional.
- For CLIP, populate **5' Barcode** when present.
- 5' barcode may be:
  - UMI-like pattern (`NNNN`, `NNNNNNNN`, etc.)
  - explicit nucleotide string (`ACGT...`)

### GEO / publication fields
- **GEO ID**: GSM accession (for example `GSM1234567`) from GEO sample records.
- **Pubmed id**: PubMed identifier from the linked publication.

### Organism/source fields
- **Organism / Species codes**:
  - `Hs` = Homo sapiens
  - `Mm` = Mus musculus
  - `Gg` = Gallus gallus
- **Source**: species + cell/tissue context as represented in the template expectations.

### CLIP target field
- **Protein (Purification Target)**: prefer the study's protein symbol used in titles/characteristics (commonly HGNC-like gene symbols or tagged constructs), not free text.
- Only use Ensembl-style identifiers when the source explicitly provides Ensembl IDs.

### RNA-Seq strandedness
- Required for RNA-Seq rows.
- Allowed values: `forward`, `reverse`, `auto`, `unstranded`.

## Normalization Guidance
- Preserve original metadata in notes if ambiguous, but normalize output values to template conventions.
- If a field is unknown, leave empty rather than guessing.
- Use consistent delimiters in sample naming (for example underscore-separated tokens).
- Ensure paired-end rows always carry both file fields.

## Validation Checklist
- [ ] Every row has `File name`.
- [ ] Paired-end rows include `File 2`; single-end rows do not.
- [ ] `Sample Name` follows required token order.
- [ ] `Type` is populated (`RNA-Seq` or `CLIP`).
- [ ] RNA-Seq rows include valid strandedness.
- [ ] CLIP rows include 5' barcode when available.
- [ ] GEO IDs use `GSM` format where provided.
- [ ] Species code is one of `Hs`, `Mm`, `Gg` when applicable.

## Output Expectations
- Produce a table compatible with `uploadsample_flowbio_v6.py` input expectations.
- Prefer explicit missing values over inferred guesses.
- If source metadata conflicts, surface a short assumptions list before finalizing.

## GEO Matrix to Run-Page Alignment
- For GEO series matrix files, treat sample-level rows as column-aligned vectors in a shared GSM order.
- Use `!Sample_geo_accession` as the canonical GSM column index map.
- The same GSM index should be used to pull aligned values from:
  - `!Sample_title` (sample title text shown on SRA run pages)
  - `!Sample_source_name_ch1`
  - `!Sample_organism_ch1`
  - `!Sample_characteristics_ch1` (multiple rows; select the row with desired label, e.g. `antibody:`)
  - `!Sample_description` (multiple rows; barcode/UMI patterns may appear in one specific description row)
  - `!Sample_instrument_model`, `!Sample_library_strategy`, and related assay fields
- When mapping SRR to GSM, prefer SRA run pages (`/sra/?term=SRR...`) because they expose both run accession and GEO accession in one record.

## Fallback for Weak SRA Titles
- If SRA run metadata returns weak or generic sample titles (or empty `LibraryName`), use the GEO sample page as primary title source:
  - `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM...`
- On GEO sample pages, parse:
  - `Title` for high-quality sample naming inputs
  - `Contact name` for `Scientist`
  - `Organization name` for `Organisation`
- This fallback is especially useful for ENCODE/ENCORE-style datasets where file-level accessions are rich but run-level title fields can be inconsistent.
- Prefer GEO `Title` over SRA `LibraryName` when the two conflict and GEO appears more descriptive.

## FLASH Project Example Conventions
Reference example files:
- `projects/flash/GSE118265_series_matrix.txt`
- `projects/flash/flashCLIP.xlsx` (and `flashCLIP.csv`)

Observed conventions to reuse when appropriate:
- Paired-end outputs may be represented as two rows with mate-specific file names (for example `_1` and `_2`) while keeping shared metadata identical.
- `Sample Name` format in this project follows a compact token style:
  - `<protein_or_target>_<species_code>_<cell_type>_<condition_or_series>_<replicate>`
- GEO IDs are GSM-level and repeated per mate when paired-end rows are split.
- CLIP rows use:
  - `Type = CLIP`
  - `Experimental Method = iCLIP`
  - `Organism` two-letter species code (`Hs` in this dataset)
  - `Cell or Tissue` as concise cell label (for example `HEK293`)
- `Purification Agent` may come from either affinity tags (`HIS-STREP`, `FLAG-STREP`) or antibody descriptors in matrix characteristics.
- 5' barcode / UMI-like pattern can be assembled from matrix tag fields and should preserve degenerate symbols (`N`, `R`, `Y`, `B`) when present.

When adopting this style in other studies:
- Keep naming token order consistent within a project.
- Prefer the project's pre-existing naming style if one already exists.
- If no style exists, fall back to the default sample-name rules in this skill.

## YEO Naming Guardrails (Training-Calibrated)

Use these rules for YEO-style large mixed cohorts where GEO/SRA metadata quality varies:

- Build sample names as:
  - `protein(if CLIP)_celltype_Hs_condition1_condition2(if present)_rep_SRR#`
- Always suffix run accession (`_SRR12345678`) for uniqueness and traceability.
- Never emit placeholder starts such as `GSM*_Hs`, `eCLIP_*`, `replicate_*`, or `UNK`.

### Source Priority for Naming Inputs
1. Curated training overrides (if provided by user for specific runs).
2. GEO sample page (`Title`, `Characteristics`, `Description`).
3. SRA run page / runinfo (`LibraryName`, `SampleName`).
4. Existing row metadata (`Notes`, `Source`, `Condition`, `Protein (Purification Target)`).

### Token Extraction Rules
- **Protein (CLIP):**
  - Prefer explicit target from `Title` / `fraction` patterns (for example `PUF60 IP`, `V5 IP`).
  - For INPUT controls use `SMInput` in metadata fields, but keep descriptive protein context in sample name when user style requires it.
  - Keep special CLIP control names such as `V5_ZFtag`.
- **Cell type / line:**
  - Normalize `HEK293T` -> `HEK293`.
  - Preserve motor neuron lineage tags (for example `CR463motor`, `CR464motor`, `Kin24motor`, `H1motor`).
  - Treat `Kin 2-4` / `Kin2ALS4` as `Kin24` family.
- **Condition(s):**
  - Use compact, informative condition tokens (`IP`, `Input`, `OE`, `Nuclear`, `Cytoplasmic`, `Insoluble`, `Recovery`, `Puromycin`, `RBDdel`).
  - For the stress/recovery study, map stress-like treatment naming to `Puromycin` to match project naming convention.
  - Keep at most two condition tokens unless user explicitly requests more.
- **Replicate:**
  - Preserve replicate token from title/description (`rep1`, `rep2`, ...).

### Metadata to Populate Alongside Sample Name
- `Type` and `Experimental Method` must be assigned before final name construction.
- Set `Source` to normalized cell context used in naming.
- Keep `Notes` with the best descriptive title text for auditability.
- Preserve `Protein (Purification Target)` and `Condition` as structured fields mirroring sample-name semantics.

### Quality Gates Before Finalizing
- Reject names that:
  - start with generic words (`eCLIP`, `replicate`, `sample`, `cell`, `type`)
  - contain placeholder fragments (`GSM`, `UNK`) when better metadata exists
  - duplicate existing names without SRR suffix
- Emit a changed-row audit table with:
  - run accession
  - old name
  - new name
  - reason/rule used

## Retrospective Process Learnings

### Execution strategy
- Process one project at a time end-to-end before mixing cohorts from other projects.
- For large projects, run in phases:
  1. matrix parsing and GSM/SRR alignment
  2. initial annotation draft
  3. focused manual training subset
  4. constrained harmonization pass
  5. final QC audit
- Never run broad name rewrites without a low-confidence filter and a change log.

### Protein name standardization
- In YEO-style CLIP datasets, protein names are usually represented as gene-symbol-like tokens in:
  - `!Sample_title`
  - `genotype:` characteristics
  - first token of curated `Notes` title segment
- Keep canonical protein tokens concise and stable (examples: `PUF60`, `ZC3HAV1`, `CNOT7`, `RNF14`, `V5ZFtag`, `L140P-PUF60-V5`).
- Avoid propagating protocol text into protein names (for example antibody catalog strings, `treatment: Bethyl ...`, `eCLIP`).
- If CLIP sample names use `RBP` placeholder and the first title token is a valid protein token, replace `RBP` with that token.

### Cell type vs protein disambiguation
- Classify as **cell type/line** if token appears in:
  - `cell line:` / `cell type:` fields
  - known cell alias set (for example `HEK293`, `H1motor`, `CR463motor`, `CR464motor`, `Kin24motor`, `MDAMB231`, `SUM149`, `MCF10A`)
- Classify as **protein target** if token appears in:
  - title position immediately preceding `replicate`
  - `genotype:` field (for engineered construct datasets)
  - `fraction: <TOKEN> IP/IN` pattern
- If a token can be both (ambiguous), resolve by precedence:
  1. explicit `cell line/cell type` assignment -> cell
  2. explicit `genotype` or `<token> replicate` title head -> protein
  3. otherwise mark for manual review, do not guess.

### Contradiction cleanup rules
- Do not force all rows to `CLIP`; preserve assay type (`CLIP` vs `RNA-Seq`) from source mapping.
- For RNA-Seq rows, do not prepend protein token in sample names.
- Never allow sample-name prefixes such as `GSM`, `eCLIP`, `replicate`, `UNK`, `rep1`, or `rep2`.
- Ensure replicate token appears at most once in each sample name.

## FlowBio `protein_target__annotation` (two underscores)

Flow stores this field on samples as **`protein_target__annotation`** (double underscore). When editing via REST, POST JSON uses the same key, for example:

```json
{"protein_target__annotation": "TARDBP-nGFP"}
```

Use this section when inferring or normalizing purification-target annotations for CLIP samples (upload tables, audits, or bulk fixes).

### Protein tags (N/C terminal, hyphen)

- **Format:** `<PROTEIN>-<tag>` with optional **N/C terminal** prefix on the tag: **`n`** or **`c`** immediately before common tags.
- **Examples:** `TARDBP-nGFP`, `TARDBP-cGFP`, `GENE-FLAG`, `GENE-V5`.
- **Hyphen** denotes a **tag fused onto** the protein target.
- **Common tags:** `FLAG`, `GFP`, `V5`, `1H4` (extend if the study uses additional standard tags).
- **Terminal marker:** use **`n`** or **`c`** when the biology is clear (N-terminal vs C-terminal fusion). If terminal is **unclear**, **omit** `n`/`c` and leave manual follow-up (example: `Srsf11-FLAG` → `SRSF11-FLAG`, not `SRSF11-nFLAG`).

### Protein mutations (before tags)

- **Format:** `<PROTEIN>:<mutation>-<tag>` when mutations are explicit.
- **Example:** `TARDBP:340del346-nGFP` means a deletion between amino acids 340 and 346, with an N-terminal GFP tag.
- **Title-style leading mutation** (example from sample titles): `L140P-PUF60-V5` normalizes to **`PUF60:L140P-V5`** (mutation scoped to the protein, then tag).

### Mutation token conventions

- **Two amino-acid letters with a number in the middle** (for example `L140P`, `R521G`) generally denotes a **single amino-acid substitution** at that residue.
- **`CTD`** denotes a **C-terminal deletion** in the protein; encode as part of the mutation token (for example `PROTEIN:CTD` or combined with other mutation syntax as the study defines).

### Special fusion naming

- Example pattern: `LIN28A_GFPNLS_Mm_nESC_Ctfusion_2` annotates as **`LIN28A-cGFP`** when the readthrough indicates **C-terminal GFP / NLS fusion** (use `cGFP` when C-terminal context is explicit).

### reCLIP (second CLIP after purification)

- If **experimental method** is **reCLIP**, use an **underscore** prefix form: **`reCLIP_<prior_target>`** (example: `reCLIP_hnRNPC`) to record that CLIP was performed **after** purification of the named target.

### Upload template vs API keys

- `uploadsample_flowbio_v6.py` maps spreadsheet column **Purification Target Annotation** to metadata key **`purification_target_annotation`** (single underscore). Flow may also expose or prefer **`protein_target__annotation`** (double underscore) in the sample JSON returned by `GET /samples/{id}`.
- Before bulk edits, confirm which key is populated in the API for your project cohort; the enrichment script `flowAPIscripts/analysis/flow_public_clip_enrich_protein_annotation.py` checks both families and records which key supplied the API value in `protein_target__annotation_api_key`.

### Workflow for bulk updates

1. **Pull** existing `protein_target__annotation` from the API (preserve non-empty values).
2. **Infer** only where the field is empty, using sample name, `purification_target`, `condition`, `notes`, and experimental method.
3. Emit a review table with columns such as: `protein_target__annotation_api`, `protein_target__annotation_inferred`, `protein_target__annotation` (merged = API if set, else inferred), and `protein_target__annotation_source`. Use `flowAPIscripts/analysis/flow_public_clip_enrich_protein_annotation.py` on `projects/flow_public_clip_samples_pull_v1.csv` (with credentials), or `--offline` for inference-only from the v1 CSV when the API is unavailable. For a **full metadata export of all public samples** (every metadata `value` plus nested `annotation` as `<key>__annotation`, plus `file_names`), run `flowAPIscripts/analysis/flow_public_samples_pull_v3.py` with credentials; output defaults to `projects/flow_public_samples_pull_v3.csv`.
4. **Do not POST** bulk edits until a human confirms merged values.

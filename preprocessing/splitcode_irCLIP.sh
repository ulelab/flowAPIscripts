#!/usr/bin/env bash
#SBATCH -p cpu
#SBATCH -t 1-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=24
#SBATCH -o logs/cutadapt_%A_%a.out
#SBATCH -e logs/cutadapt_%A_%a.err

set -euo pipefail

# --- paths ---
CFG="../config.txt"          # splitcode config (has @extract for RX etc.)
IN_DIR="./"                  # directory containing SRR*.trimmed.fq.gz
OUT_DIR="../cleaned"         # final FASTQs after trimming last 3 nt
LOG_DIR="../logs"
TMP_DIR="../tmp_splitcode"
mkdir -p "$OUT_DIR" "$LOG_DIR" "$TMP_DIR"

ASSIGN_LOG="$LOG_DIR/splitcode_assigned.csv"
: > "$ASSIGN_LOG"
echo "sample,outfile,assigned_percent" >> "$ASSIGN_LOG"

prepend_outfile_to_mapping() {
  local outfile="$1" mapping_in="$2" mapping_out="$3"
  { printf "outfile: %s\n" "$outfile"; cat "$mapping_in"; } > "$mapping_out"
}

# NOTE: match ONLY adapter-trimmed inputs; avoid any *untrimmed*.fq.gz files.
for gz in "$IN_DIR"/SRR*.trimmed.fq.gz; do
  [[ -e "$gz" ]] || { echo "No SRR*.trimmed.fq.gz files found in $IN_DIR"; break; }

  base="$(basename "$gz" .trimmed.fq.gz)"              # e.g., SRR29661270
  tmp_fastq="$TMP_DIR/${base}.splitcode.fastq"         # pass 1 output (plain FASTQ)
  out_fastq="$OUT_DIR/${base}_cleaned.fastq"           # pass 2 output (plain; gzip after)
  out_fastq_gz="${out_fastq}.gz"
  tmp_stderr="$(mktemp)"
  tmp_map="$(mktemp)"
  per_file_map="${OUT_DIR}/${base}_mapping.txt"

  echo "[INFO] Processing $gz → $out_fastq_gz (splitcode extract, then trim last 3 nt)"

  # PASS 1: splitcode extraction only (no trimming).
  splitcode \
    -c "$CFG" \
    -N 1 \
    --x-names \
    --assign \
    --mapping="$tmp_map" \
    --no-outb \
    -o "$tmp_fastq" \
    "$gz" \
    2> "$tmp_stderr"

  # PASS 2: trim last 3 nt with splitcode CLI (equivalent to header '@trim-3 3')
  splitcode \
    --trim-3 3 \
    -N 1 \
    -o "$out_fastq" \
    --keep-com \
    "$tmp_fastq" \
    1> /dev/null

  # Append Assigned% from PASS 1
  assigned=$(grep -i -E 'assigned' "$tmp_stderr" | grep -Eo '[0-9]+(\.[0-9]+)?%' || true)
  [[ -n "${assigned:-}" ]] || assigned="NA"
  echo "${base},${base}_cleaned.fastq.gz,${assigned}" >> "$ASSIGN_LOG"

  # Save annotated mapping with final outfile name
  prepend_outfile_to_mapping "${base}_cleaned.fastq.gz" "$tmp_map" "$per_file_map"

  # Compress final FASTQ; use Slurm CPU allocation if present
  pigz -p "${SLURM_CPUS_PER_TASK:-4}" -f "$out_fastq"

  rm -f "$tmp_stderr" "$tmp_map" "$tmp_fastq"
done

echo "[INFO] Done. Assigned% summary at: $ASSIGN_LOG"

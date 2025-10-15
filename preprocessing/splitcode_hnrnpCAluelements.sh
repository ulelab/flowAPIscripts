#!/usr/bin/env bash
set -euo pipefail

CFG="./config.txt"
IN_DIR="./"
OUT_DIR="trimmings"
LOG_DIR="logs"
mkdir -p "$OUT_DIR" "$LOG_DIR"

ASSIGN_LOG="$LOG_DIR/splitcode_assigned.csv"
: > "$ASSIGN_LOG"
echo "sample,outfile,assigned_percent,umi_extracted_percent,n_processed,n_assigned,n_umi" >> "$ASSIGN_LOG"

keep_rx_reads() {
  local in_fastq="$1" out_fastq="$2"
  # keep 4-line blocks whose header contains RX:Z:
  awk 'NR%4==1 {keep = ($0 ~ /RX:Z:/)} { if (keep) print }' "$in_fastq" > "$out_fastq"
}

for gz in "$IN_DIR"/ERR*.fastq.gz; do
  [[ -e "$gz" ]] || { echo "No ERR*.fastq.gz files found in $IN_DIR"; break; }

  base="$(basename "$gz" .fastq.gz)"
  out_fastq="${OUT_DIR}/${base}_trimmed.fastq"
  tmp_map="$(mktemp)"
  tmp_summary="$(mktemp)"
  per_file_map="${OUT_DIR}/${base}_mapping.txt"

  echo "[INFO] Processing $gz → $out_fastq"

  splitcode \
    -c "$CFG" \
    -N 1 \
    --x-names \
    --assign \
    --mapping="$tmp_map" \
    --summary="$tmp_summary" \
    --no-outb \
    -o "$out_fastq" \
    "$gz"

  # Parse JSON summary from stdin to avoid argv/heredoc pitfalls
  if [[ -s "$tmp_summary" ]]; then
    read -r n_processed n_assigned n_umi <<<"$(
      python3 - <<'PY' < "$tmp_summary"
import sys, json
try:
    j=json.load(sys.stdin)
except Exception:
    print("0 0 0"); sys.exit(0)
n_proc=int(j.get("n_processed",0))
n_ass =int(j.get("n_assigned",0))
n_umi =0
for x in j.get("extraction_info",[]):
    if x.get("name")=="umi":
        n_umi=int(x.get("n_reads",0)); break
print(n_proc, n_ass, n_umi)
PY
    )"
  else
    n_processed=0; n_assigned=0; n_umi=0
  fi

  # Percentages
  if [[ "${n_processed:-0}" -gt 0 ]]; then
    assigned_pct=$(awk -v a="$n_assigned" -v p="$n_processed" 'BEGIN{printf "%.1f%%", (100.0*a/p)}')
    umi_pct=$(awk -v u="$n_umi" -v p="$n_processed" 'BEGIN{printf "%.1f%%", (100.0*u/p)}')
  else
    assigned_pct="NA"; umi_pct="NA"
  fi
  echo "${base},${base}_trimmed.fastq,${assigned_pct},${umi_pct},${n_processed},${n_assigned},${n_umi}" >> "$ASSIGN_LOG"

  # Save mapping with the outfile name at top
  { printf "outfile: %s\n" "${base}_trimmed.fastq"; cat "$tmp_map"; } > "$per_file_map"

  # OPTIONAL: produce a filtered file that drops reads without an RX tag in the header
  filtered_fastq="${OUT_DIR}/${base}_trimmed.umiOnly.fastq"
  keep_rx_reads "$out_fastq" "$filtered_fastq"

  nohup pigz -p "$(nproc)" -f "$out_fastq" "$filtered_fastq" >/dev/null 2>&1 &

  rm -f "$tmp_map" "$tmp_summary"
done

echo "[INFO] Done launching pigz jobs. Active compressions:"
pgrep -a pigz || true

echo "[INFO] Assigned% summary at: $ASSIGN_LOG"

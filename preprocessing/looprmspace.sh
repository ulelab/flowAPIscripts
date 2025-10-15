#!/usr/bin/env bash
#SBATCH -p cpu
#SBATCH -t 1-00:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=12
#SBATCH --array=1-116
#SBATCH -o logs/rmspace_%A_%a.out
#SBATCH -e logs/rmspace_%A_%a.err

set -euo pipefail
shopt -s nullglob

IN_DIR="${IN_DIR:-.}"
SCRIPT="${SCRIPT:-../removespace.py}"
SAMPLESHEET="${SAMPLESHEET:-../irCLIP.tsv}"
mkdir -p logs

if [[ ! -f "$SAMPLESHEET" ]]; then
  echo "[ERROR] Samplesheet not found: $SAMPLESHEET" >&2
  exit 1
fi

# Read first column from samplesheet, ignore blank lines and lines starting with '#'.
# Resolve bare basenames relative to IN_DIR. Only keep files that exist.
mapfile -t FILES < <(awk -F"\t" 'NR>1 && NF && $1 !~ /^#/ {print $1}' "$SAMPLESHEET" | sed 's/\r$//' | while IFS= read -r entry; do
  [[ -z "$entry" ]] && continue
  if [[ "$entry" == /* || "$entry" == .* || "$entry" == */* ]]; then
    path="$entry"
  else
    path="$IN_DIR/$entry"
  fi
  if [[ -f "$path" ]]; then
    printf '%s\n' "$path"
  else
    printf '[WARN] Missing file: %s\n' "$path" >&2
  fi
done)

N=${#FILES[@]}
if (( N == 0 )); then
  echo "[ERROR] No valid files found from first column of $SAMPLESHEET" >&2
  exit 1
fi

# Convert SLURM array id (1-based) to index in FILES (0-based)
IDX=$(( SLURM_ARRAY_TASK_ID - 1 ))
if (( IDX < 0 )); then
  echo "[ERROR] SLURM_ARRAY_TASK_ID must be >= 1" >&2
  exit 1
fi
if (( IDX >= N )); then
  echo "[INFO] Array index $SLURM_ARRAY_TASK_ID exceeds valid file count $N; nothing to do." >&2
  exit 0
fi

IN="${FILES[$IDX]}"
BASE="$(basename "$IN")"

echo "[INFO] ($((IDX+1))/$N) Processing $BASE"
echo "[INFO] Using script: $SCRIPT"

# Call removespace.py, which handles its own output naming
python3 "$SCRIPT" "$IN"

echo "[OK] Finished processing: $BASE"

#!/usr/bin/env bash
#SBATCH -p cpu
#SBATCH -t 1-00:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-246
#SBATCH -o logs/cutadapt_%A_%a.out
#SBATCH -e logs/cutadapt_%A_%a.err

set -euo pipefail

# --- paths ---
INPUT_DIR="${INPUT_DIR:-$(pwd)/fastq_files}"
OUT_DIR="${OUT_DIR:-$(pwd)/trimmed}"
LOG_DIR="${LOG_DIR:-$(pwd)/logs}"
mkdir -p "$OUT_DIR" "$LOG_DIR"

# --- conda env: activate 'cutadapt' ---

source /scratch/prj/ppn_rnp_networks/users/mike.jones/software/mambaforge/etc/profile.d/conda.sh
conda activate /scratch/prj/ppn_rnp_networks/users/mike.jones/software/mambaforge/envs/cutadapt

# --- pick the SRR file for this array index ---
# Accept .fastq.gz or .fq.gz
mapfile -t FILES < <(find "$INPUT_DIR" -maxdepth 1 -type f \( -name 'SRR*.fastq.gz' -o -name 'SRR*.fq.gz' \) | sort)
N=${#FILES[@]}
IDX=$((SLURM_ARRAY_TASK_ID - 1))

if (( N == 0 )); then
  echo "No SRR*.fastq.gz or SRR*.fq.gz files found in $INPUT_DIR" >&2
  exit 1
fi
if (( IDX >= N )); then
  echo "Array index $SLURM_ARRAY_TASK_ID has no corresponding file (only $N files)." >&2
  exit 0
fi

IN="${FILES[$IDX]}"
BASE="$(basename "$IN")"
SAMPLE="${BASE%.fastq.gz}"
SAMPLE="${SAMPLE%.fq.gz}"

UNTRIMMED="$OUT_DIR/${SAMPLE}.untrimmed.fq.gz"
TRIMMED="$OUT_DIR/${SAMPLE}.trimmed.fq.gz"
REPORT="$LOG_DIR/${SAMPLE}.trimmed.report.txt"

echo "[$(date)] Processing: $IN"
echo "Outputs: $TRIMMED , $UNTRIMMED"
echo "Report:  $REPORT"

# --- run cutadapt (your settings) ---
cutadapt \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
  --overlap 5 \
  -e 0.2 \
  --minimum-length 46 \
  --untrimmed-output "$UNTRIMMED" \
  -o "$TRIMMED" \
  "$IN" > "$REPORT" 2>&1

echo "[$(date)] Done."
#!/bin/bash
#SBATCH --job-name=uploadsample
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-29
#SBATCH --output=/scratch/prj/ppn_rnp_networks/users/mike.jones/logs/uploadsample_simple_%A_%a.out
#SBATCH --error=/scratch/prj/ppn_rnp_networks/users/mike.jones/logs/uploadsample_simple_%A_%a.err

set -euo pipefail

ENV_PREFIX="/scratch/prj/ppn_rnp_networks/users/mike.jones/software/mambaforge"
source "${ENV_PREFIX}/etc/profile.d/conda.sh"
conda activate flowbio

XLSX_FILE="/scratch/prj/ppn_rnp_networks/users/mike.jones/data/speckles/SraRunTable_from_template.xlsx"

ROW_SPEC=$(python3 - "${XLSX_FILE}" "${SLURM_ARRAY_TASK_ID}" 1 <<'PY'
import sys
import pandas as pd

xlsx_path = sys.argv[1]
task_id = int(sys.argv[2])
samples_per_task = int(sys.argv[3])

# Calculate the starting row index for this task (1-based)
# One sample per task: Task 1 -> row 1, Task 2 -> row 2, ...
start_idx = (task_id - 1) * samples_per_task + 1
end_idx = start_idx + samples_per_task - 1

# Read Excel file
df = pd.read_excel(xlsx_path)
num_rows = len(df)

rows = []
for idx in range(1, num_rows + 1):
    if start_idx <= idx <= end_idx:
        rows.append(str(idx))

if rows:
    print(",".join(rows))
PY
)

if [[ -z "${ROW_SPEC}" ]]; then
  echo "No rows assigned to array task ${SLURM_ARRAY_TASK_ID}; exiting."
  exit 0
fi

python3 uploadsample_flowbio_v5.py \
  -p 392381750398848919 \
  -r "${ROW_SPEC}" \
  "${XLSX_FILE}"

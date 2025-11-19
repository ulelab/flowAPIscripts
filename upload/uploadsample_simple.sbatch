#!/bin/bash
#SBATCH --job-name=uploadsample
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-50
#SBATCH --output=/scratch/prj/ppn_rnp_networks/users/mike.jones/logs/uploadsample_simple_%A_%a.out
#SBATCH --error=/scratch/prj/ppn_rnp_networks/users/mike.jones/logs/uploadsample_simple_%A_%a.err

set -euo pipefail

ENV_PREFIX="/scratch/prj/ppn_rnp_networks/users/mike.jones/software/mambaforge"
source "${ENV_PREFIX}/etc/profile.d/conda.sh"
conda activate flowbio

TSV_FILE="/scratch/prj/ppn_rnp_networks/users/mike.jones/data/TDP43QuantSeq/TDPcondensateRNA.tsv"
BASE_DIR="/scratch/prj/ppn_rnp_networks/users/mike.jones/data/TDP43QuantSeq"

ROW_SPEC=$(python3 - "${TSV_FILE}" "${SLURM_ARRAY_TASK_ID}" 50 <<'PY'
import csv
import sys

tsv_path = sys.argv[1]
task_id = int(sys.argv[2])
stride = int(sys.argv[3])

rows = []
with open(tsv_path, newline="") as handle:
    reader = csv.reader(handle, delimiter="\t")
    header = next(reader, None)
    for idx, _ in enumerate(reader, start=1):
        if idx < task_id:
            continue
        if (idx - task_id) % stride == 0:
            rows.append(str(idx))

if rows:
    print(",".join(rows))
PY
)

if [[ -z "${ROW_SPEC}" ]]; then
  echo "No rows assigned to array task ${SLURM_ARRAY_TASK_ID}; exiting."
  exit 0
fi

python3 ../../scripts/uploadsample_static.py \
  "${TSV_FILE}" \
  --rows "${ROW_SPEC}" \
  --base-dir "${BASE_DIR}"


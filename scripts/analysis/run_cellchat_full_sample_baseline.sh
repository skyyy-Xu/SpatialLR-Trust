#!/bin/bash
#SBATCH -J slr_cellchat_full1
#SBATCH -o runs/%x.%j.out
#SBATCH -e runs/%x.%j.err
#SBATCH -p normal
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH -t 04:00:00

set -euo pipefail
umask 027

PROJECT_DIR="${PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
RUN_ID=${RUN_ID:-20260712_1701_cellchat-full-sample-calibration}
SAMPLE_ID=${SAMPLE_ID:-GSM9060732_AdjIII-0019}
OUTPUT_PREFIX=${OUTPUT_PREFIX:-cellchat_full_gsm9060732}
ENV_DIR="${PROJECT_DIR}/envs/spatiallr-r-cellchat-cytosignal"
R_SITE_LIBRARY="${ENV_DIR}/lib/R/site-library"
RSCRIPT="${ENV_DIR}/bin/Rscript"
MARKER="${ENV_DIR}/.spatiallr_r_env_complete"
RUN_DIR="${PROJECT_DIR}/runs/${RUN_ID}"
OUT_DIR="${PROJECT_DIR}/results/task_c_cellchat_full_sample_baseline"

cd "${PROJECT_DIR}"
mkdir -p "${RUN_DIR}" "${OUT_DIR}"
chmod 750 "${RUN_DIR}" "${OUT_DIR}"
export RUN_ID R_LIBS_USER="${R_SITE_LIBRARY}"
export PATH="${ENV_DIR}/bin:${PATH}"

{
  echo "RUN_ID=${RUN_ID}"
  echo "SLURM_JOB_ID=${SLURM_JOB_ID:-manual}"
  echo "HOST=$(hostname)"
  echo "START=$(date -Is)"
  echo "PROJECT_DIR=${PROJECT_DIR}"
  echo "ENV_DIR=${ENV_DIR}"
  echo "RSCRIPT=${RSCRIPT}"
  echo "SAMPLE_ID=${SAMPLE_ID}"
  echo "OUTPUT_PREFIX=${OUTPUT_PREFIX}"
  echo "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-1}"
  echo "SLURM_MEM_PER_NODE=${SLURM_MEM_PER_NODE:-unknown}"
} | tee "${RUN_DIR}/environment.txt"

if [[ ! -x "${RSCRIPT}" || ! -d "${R_SITE_LIBRARY}" || ! -e "${MARKER}" ]]; then
  echo "ERROR: completed project-local CellChat environment is unavailable." >&2
  exit 2
fi

"${RSCRIPT}" scripts/server/run_cellchat_full_sample_baseline.R \
  --project_dir "${PROJECT_DIR}" \
  --sample_id "${SAMPLE_ID}" \
  --output_prefix "${OUTPUT_PREFIX}" \
  --min_cells "${MIN_CELLS:-10}" \
  --trim "${TRIM:-0.1}" \
  --seed "${SEED:-20260712}"

cat > "${RUN_DIR}/outputs_manifest.tsv" <<MANIFEST_EOF
path	type	description
results/task_c_cellchat_full_sample_baseline/${OUTPUT_PREFIX}_raw_subsetCommunication.tsv	text	Raw CellChat full-expression interactions
results/task_c_cellchat_full_sample_baseline/${OUTPUT_PREFIX}_candidates.tsv	text	Standardized CellChat full-expression candidates
results/task_c_cellchat_full_sample_baseline/${OUTPUT_PREFIX}_summary.json	json	Full-expression resource calibration summary
runs/${RUN_ID}/environment.txt	text	Slurm execution metadata
runs/slr_cellchat_full1.${SLURM_JOB_ID:-manual}.out	log	Slurm standard output
runs/slr_cellchat_full1.${SLURM_JOB_ID:-manual}.err	log	Slurm standard error
MANIFEST_EOF
chmod -R o-rwx "${RUN_DIR}" "${OUT_DIR}"
echo "END=$(date -Is)"
echo "DONE"

#!/bin/bash
#SBATCH --job-name=conn-cmp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=48:00:00
#SBATCH --array=0-6
## Logs (create directory before first submit: mkdir -p Connectivity/slurm/logs)
#SBATCH --output=Connectivity/slurm/logs/%x_%A_%a.out
#SBATCH --error=Connectivity/slurm/logs/%x_%A_%a.err

# Run single-vs-bootstrap comparison for one coral family per array task.
# Families match: connectivity_acroporidae.nc, ..._agariciidae.nc, etc.
#
# Submit from the repository ROOT (parent of Connectivity/), e.g.:
#   export CONNECTIVITY_MATRICES_DIR=/path/to/dir/with/netcdfs
#   mkdir -p Connectivity/slurm/logs
#   sbatch Connectivity/slurm/run_family_array.slurm
#
# Optional:
#   export CONNECTIVITY_REPO_ROOT=/path/to/Connectivity_analysis   # if not using submit dir
#   sbatch --export=ALL,CONNECTIVITY_MATRICES_DIR=/data/nc Connectivity/slurm/run_family_array.slurm

set -euo pipefail

if [[ -z "${CONNECTIVITY_MATRICES_DIR:-}" ]]; then
  echo "ERROR: Set CONNECTIVITY_MATRICES_DIR to the directory containing connectivity_*.nc files."
  exit 1
fi

MATRIX_DIR="$(readlink -f "${CONNECTIVITY_MATRICES_DIR}")"
REPO_ROOT="${CONNECTIVITY_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
REPO_ROOT="$(readlink -f "${REPO_ROOT}")"
CONN_DIR="${REPO_ROOT}/Connectivity"

if [[ ! -d "${CONN_DIR}" ]]; then
  echo "ERROR: Connectivity code directory not found: ${CONN_DIR}"
  exit 1
fi

# Order must match your file set (bootstrap connectivity_<family>.nc)
FAMILIES=(
  acroporidae
  agariciidae
  diploastraeidae
  euphylliidae
  lobophyllidae
  merulinidae
  poritidae
)

FAMILY="${FAMILIES[${SLURM_ARRAY_TASK_ID}]}"
BOOT="${MATRIX_DIR}/connectivity_${FAMILY}.nc"
SINGLE="${MATRIX_DIR}/connectivity_${FAMILY}_single.nc"
OUT="${CONN_DIR}/outputs_${FAMILY}"

echo "=========================================="
echo "Task ID     : ${SLURM_ARRAY_TASK_ID}"
echo "Family      : ${FAMILY}"
echo "Bootstrap   : ${BOOT}"
echo "Single      : ${SINGLE}"
echo "Output dir  : ${OUT}"
echo "Host        : $(hostname)"
echo "Start       : $(date -Is)"
echo "=========================================="

if [[ ! -f "${BOOT}" ]]; then
  echo "ERROR: Bootstrap file missing: ${BOOT}"
  exit 2
fi
if [[ ! -f "${SINGLE}" ]]; then
  echo "ERROR: Single-matrix file missing: ${SINGLE}"
  exit 2
fi

# --- Optional: HPC modules (uncomment and edit for your site) ---
# module purge
# module load miniforge3/24.3  # example
# source "$(conda info --base)/etc/profile.d/conda.sh"
# conda activate gbr-connectivity

cd "${CONN_DIR}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

python run_analysis.py \
  --bootstrap-path "${BOOT}" \
  --single-path "${SINGLE}" \
  --output-dir "${OUT}"

echo "Finished: $(date -Is)"

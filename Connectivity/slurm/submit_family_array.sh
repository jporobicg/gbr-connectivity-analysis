#!/usr/bin/env bash
# Helper: submit the family array job from the repository root.
# Usage:
#   ./Connectivity/slurm/submit_family_array.sh /path/to/directory/with/connectivity_*.nc
#
# Or from repo root:
#   bash Connectivity/slurm/submit_family_array.sh "$PWD/datasets/connectivity_matrices"

set -euo pipefail

MATRIX_DIR="${1:?Pass the directory containing connectivity_<family>.nc and connectivity_<family>_single.nc files}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

export CONNECTIVITY_MATRICES_DIR="$(readlink -f "${MATRIX_DIR}")"
export CONNECTIVITY_REPO_ROOT="${REPO_ROOT}"

mkdir -p "${REPO_ROOT}/Connectivity/slurm/logs"

cd "${REPO_ROOT}"
sbatch "${REPO_ROOT}/Connectivity/slurm/run_family_array.slurm"

echo "Submitted from ${REPO_ROOT} with CONNECTIVITY_MATRICES_DIR=${CONNECTIVITY_MATRICES_DIR}"

#!/bin/bash
#SBATCH --job-name=hawkes_bivar_est
#SBATCH --array=1-400
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=24:00:00

set -euo pipefail

CODE_DIR="${HAWKES_BIVAR_CODE_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
RUN_DIR="${HAWKES_BIVAR_RUN_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
CODE_DIR="$(cd "$CODE_DIR" && pwd)"
mkdir -p "$RUN_DIR/logs" "$RUN_DIR/bivar_results" "$RUN_DIR/tmp"
export TMPDIR="${TMPDIR:-$RUN_DIR/tmp}"

# Uncomment/adapt for your cluster if needed.
module load R

ARRAY_COUNT="${ARRAY_COUNT:-${SLURM_ARRAY_TASK_COUNT:-400}}"
ARRAY_ID="${SLURM_ARRAY_TASK_ID:-1}"

Rscript "$CODE_DIR/run_bivar_hpc.R" \
  --base-dir "$CODE_DIR" \
  --array-id "$ARRAY_ID" \
  --array-count "$ARRAY_COUNT" \
  --workers "${SLURM_CPUS_PER_TASK:-1}" \
  --outdir "$RUN_DIR/bivar_results" \
  --T-grid 500,1000,2000,4000,8000,16000 \
  --reps-per-T 2000 \
  --nstart-mle 20 \
  --nstart-gmm 30 \
  --nstart-aug 30 \
  --run-aug TRUE \
  --aug-s 0.4,0.4 \
  --aug-degree 1 \
  --aug-quad-hmax 0.5 \
  --aug-ridge-rel 1e-8 \
  --compute-info FALSE \
  --seed 20260525

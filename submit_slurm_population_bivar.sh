#!/bin/bash
#SBATCH --job-name=hawkes_bivar_pop
#SBATCH --array=1-64
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00

set -euo pipefail

CODE_DIR="${HAWKES_BIVAR_CODE_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
RUN_DIR="${HAWKES_BIVAR_RUN_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
CODE_DIR="$(cd "$CODE_DIR" && pwd)"
mkdir -p "$RUN_DIR/logs" "$RUN_DIR/bivar_population" "$RUN_DIR/tmp"
export TMPDIR="${TMPDIR:-$RUN_DIR/tmp}"

# Uncomment/adapt for your cluster if needed.
module load R

ARRAY_COUNT="${ARRAY_COUNT:-${SLURM_ARRAY_TASK_COUNT:-64}}"
ARRAY_ID="${SLURM_ARRAY_TASK_ID:-1}"

Rscript "$CODE_DIR/run_population_godambe_bivar.R" \
  --base-dir "$CODE_DIR" \
  --array-id "$ARRAY_ID" \
  --array-count "$ARRAY_COUNT" \
  --workers "${SLURM_CPUS_PER_TASK:-1}" \
  --outdir "$RUN_DIR/bivar_population" \
  --B-pop 128 \
  --T-pop 250000 \
  --n-eval 250000 \
  --aug-s 0.4,0.4 \
  --aug-degree 1 \
  --seed 20260525

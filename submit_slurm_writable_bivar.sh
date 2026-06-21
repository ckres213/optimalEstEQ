#!/bin/bash
# Submit the full bivariate study from a writable run directory.
# Usage:
#   ./submit_slurm_writable_bivar.sh /path/to/writable/run_dir
# If no run_dir is supplied, defaults to $SCRATCH/hawkes_bivar_godambe_<timestamp>,
# then $HOME/hawkes_bivar_godambe_<timestamp> if SCRATCH is unset.

set -euo pipefail

CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAMP="$(date +%Y%m%d_%H%M%S)"
DEFAULT_PARENT="${SCRATCH:-$HOME}"
RUN_DIR="${1:-${HAWKES_BIVAR_RUN_DIR:-$DEFAULT_PARENT/hawkes_bivar_godambe_$STAMP}}"
RUN_DIR="$(mkdir -p "$RUN_DIR" && cd "$RUN_DIR" && pwd)"
mkdir -p "$RUN_DIR/logs" "$RUN_DIR/tmp" "$RUN_DIR/bivar_population" "$RUN_DIR/bivar_results" "$RUN_DIR/bivar_figures" "$RUN_DIR/bivar_figures_paper"

export HAWKES_BIVAR_CODE_DIR="$CODE_DIR"
export HAWKES_BIVAR_RUN_DIR="$RUN_DIR"

common=(--chdir="$RUN_DIR" --export=ALL,HAWKES_BIVAR_CODE_DIR="$CODE_DIR",HAWKES_BIVAR_RUN_DIR="$RUN_DIR")

pop_jid=$(sbatch --parsable "${common[@]}" \
  --output="$RUN_DIR/logs/pop_%A_%a.out" \
  --error="$RUN_DIR/logs/pop_%A_%a.err" \
  "$CODE_DIR/submit_slurm_population_bivar.sh")

est_jid=$(sbatch --parsable "${common[@]}" \
  --output="$RUN_DIR/logs/est_%A_%a.out" \
  --error="$RUN_DIR/logs/est_%A_%a.err" \
  "$CODE_DIR/submit_slurm_estimators_bivar.sh")

comb_jid=$(sbatch --parsable "${common[@]}" \
  --dependency="afterok:${pop_jid}:${est_jid}" \
  --output="$RUN_DIR/logs/combine_%j.out" \
  --error="$RUN_DIR/logs/combine_%j.err" \
  "$CODE_DIR/submit_slurm_combine_figures_bivar.sh")

cat <<MSG
Submitted bivariate Hawkes Godambe study.
  code dir: $CODE_DIR
  run dir:  $RUN_DIR
  population job: $pop_jid
  estimator job:  $est_jid
  combine job:    $comb_jid  (after both arrays finish successfully)

After completion, check:
  $CODE_DIR/check_bivar_outputs.sh $RUN_DIR/bivar_results $RUN_DIR/bivar_population $RUN_DIR/bivar_figures $RUN_DIR/bivar_figures_paper
MSG

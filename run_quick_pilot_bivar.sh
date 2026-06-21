#!/bin/bash
set -euo pipefail

CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_DIR="${HAWKES_BIVAR_RUN_DIR:-$CODE_DIR}"
mkdir -p "$RUN_DIR/pilot_population" "$RUN_DIR/pilot_results" "$RUN_DIR/pilot_figures" "$RUN_DIR/pilot_figures_paper" "$RUN_DIR/tmp"
export TMPDIR="${TMPDIR:-$RUN_DIR/tmp}"

Rscript "$CODE_DIR/run_population_godambe_bivar.R" \
  --base-dir "$CODE_DIR" \
  --array-id 1 --array-count 1 --workers 1 \
  --outdir "$RUN_DIR/pilot_population" \
  --B-pop 2 --T-pop 5000 --n-eval 5000 \
  --aug-s 0.4,0.4 --aug-degree 1 \
  --seed 20260525 --overwrite TRUE

Rscript "$CODE_DIR/run_bivar_hpc.R" \
  --base-dir "$CODE_DIR" \
  --array-id 1 --array-count 1 --workers 1 \
  --outdir "$RUN_DIR/pilot_results" \
  --T-grid 300,600 --reps-per-T 4 \
  --nstart-mle 4 --nstart-gmm 4 --nstart-aug 4 \
  --run-aug TRUE --aug-degree 1 --aug-quad-hmax 0.5 --aug-s 0.4,0.4 \
  --aug-ridge-rel 1e-8 \
  --compute-info FALSE \
  --seed 20260525 --overwrite TRUE

Rscript "$CODE_DIR/combine_bivar_results.R" \
  --base-dir "$CODE_DIR" \
  --indir "$RUN_DIR/pilot_results" \
  --popdir "$RUN_DIR/pilot_population" \
  --outdir "$RUN_DIR/pilot_results" \
  --B-boot 50

Rscript "$CODE_DIR/make_bivar_figures.R" \
  --indir "$RUN_DIR/pilot_results" \
  --outdir "$RUN_DIR/pilot_figures"

Rscript "$CODE_DIR/paper_bivar_four_plots.R" \
  --indir "$RUN_DIR/pilot_results" \
  --popdir "$RUN_DIR/pilot_population" \
  --outdir "$RUN_DIR/pilot_figures_paper" \
  --prefix paper_bivar_pilot \
  --tmin 300


Rscript "$CODE_DIR/make_article_bivar_pair.R" \
  --indir "$RUN_DIR/pilot_results" \
  --popdir "$RUN_DIR/pilot_population" \
  --outdir "$RUN_DIR/pilot_figures_paper" \
  --prefix paper_bivar_pilot \
  --tmin 300

"$CODE_DIR/check_bivar_outputs.sh" "$RUN_DIR/pilot_results" "$RUN_DIR/pilot_population" "$RUN_DIR/pilot_figures" "$RUN_DIR/pilot_figures_paper"

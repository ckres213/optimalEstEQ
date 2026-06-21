#!/bin/bash
#SBATCH --job-name=hawkes_bivar_combine
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00

set -euo pipefail

CODE_DIR="${HAWKES_BIVAR_CODE_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
RUN_DIR="${HAWKES_BIVAR_RUN_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
CODE_DIR="$(cd "$CODE_DIR" && pwd)"
mkdir -p "$RUN_DIR/logs" "$RUN_DIR/bivar_results" "$RUN_DIR/bivar_figures" "$RUN_DIR/bivar_figures_paper" "$RUN_DIR/tmp"
export TMPDIR="${TMPDIR:-$RUN_DIR/tmp}"

# Adapt this line to your cluster if needed, e.g. module load R/4.3.2
module load R

Rscript "$CODE_DIR/combine_bivar_results.R" \
  --base-dir "$CODE_DIR" \
  --indir "$RUN_DIR/bivar_results" \
  --popdir "$RUN_DIR/bivar_population" \
  --outdir "$RUN_DIR/bivar_results" \
  --B-boot 500

Rscript "$CODE_DIR/make_bivar_figures.R" \
  --indir "$RUN_DIR/bivar_results" \
  --outdir "$RUN_DIR/bivar_figures"

# The original paper-facing set retained for diagnostics/appendix material.
Rscript "$CODE_DIR/paper_bivar_four_plots.R" \
  --indir "$RUN_DIR/bivar_results" \
  --popdir "$RUN_DIR/bivar_population" \
  --outdir "$RUN_DIR/bivar_figures_paper" \
  --prefix paper_bivar \
  --tmin 1000

# The two manuscript figures referenced in main.tex, plus CI/eigenvalue diagnostics.
Rscript "$CODE_DIR/make_final_bivar_paper_figures.R" \
  --indir "$RUN_DIR/bivar_results" \
  --popdir "$RUN_DIR/bivar_population" \
  --outdir "$RUN_DIR/bivar_figures_paper" \
  --prefix paper_bivar \
  --tmin 1000

# Full-width side-by-side article figure with fonts sized for width=\textwidth.
Rscript "$CODE_DIR/make_article_bivar_pair.R" \
  --indir "$RUN_DIR/bivar_results" \
  --popdir "$RUN_DIR/bivar_population" \
  --outdir "$RUN_DIR/bivar_figures_paper" \
  --prefix paper_bivar \
  --tmin 1000

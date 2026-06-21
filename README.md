This repository contains the code for the numerical experiment described in Section 5 of *Optimal Estimating Equations for Compact-Memory Hawkes Processes* by Louis Davis and Conor Kresin.

# Package for numerical experiment

Reproducibility code for a bivariate linear Hawkes simulation study comparing maximum likelihood with just-identified and over-identified GMM estimators. The repository contains the R/Rcpp implementation, Slurm array scripts, result aggregation, diagnostics, and manuscript figure generation.

## Requirements

- R with the `Rcpp` package
- A C++ compiler compatible with Rcpp
- Slurm (`sbatch`) for the full study

## Quick smoke test

```bash
chmod +x *.sh *.py
module load R
Rscript -e 'if (!requireNamespace("Rcpp", quietly=TRUE)) install.packages("Rcpp", repos="https://cloud.r-project.org")'
./run_quick_pilot_bivar.sh

## Key files

```text
hawkes_bivar_fast.cpp                  Rcpp implementation
bivar_core.R                           shared R utilities and estimator routines
run_population_godambe_bivar.R         population Godambe array runner
run_bivar_hpc.R                        Monte Carlo estimator array runner
combine_bivar_results.R                combines array CSVs and creates diagnostics CSVs
make_bivar_figures.R                   standard PDF diagnostics
make_final_bivar_paper_figures.R       original paper-facing individual figures
paper_bivar_four_plots.R               additional paper/appendix diagnostics
make_article_bivar_pair.R              final full-width side-by-side manuscript figure
submit_slurm_population_bivar.sh       Slurm population array script
submit_slurm_estimators_bivar.sh       Slurm estimator array script
submit_slurm_combine_figures_bivar.sh  Slurm combine-and-plot script
submit_slurm_writable_bivar.sh         wrapper for a writable run directory
run_quick_pilot_bivar.sh               small smoke test
check_bivar_outputs.sh                 output checker
patch_tex_bivar_side_by_side.py        patches main TeX to use the combined figure
focused_bivariate_simulation_section.tex manuscript section snippet
```

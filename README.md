Bivariate Hawkes Godambe simulation

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

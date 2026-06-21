# HawkesAnnalStudy clean rerun package

This is the cleaned minimal rerun package for the bivariate Hawkes Godambe simulation study. It contains source code and Slurm scripts only; old results, logs, backup files, and nested duplicate folders are intentionally excluded.

The package generates the manuscript-ready full-width two-panel figure:

```text
bivar_figures_paper/paper_bivar_efficiency_pair_Tge1000.pdf
bivar_figures_paper/paper_bivar_efficiency_pair_Tge1000.png
```

The two-panel figure uses larger article-readable fonts, keeps the panels wide, and puts the legends in a dedicated bottom strip instead of inside the plotting regions.

## Main run commands

```bash
unzip HawkesAnnalStudy_clean_final.zip
cd HawkesAnnalStudy_clean_final
chmod +x *.sh *.R *.py

module load R
Rscript -e 'if (!requireNamespace("Rcpp", quietly=TRUE)) install.packages("Rcpp", repos="https://cloud.r-project.org")'

export HAWKES_BIVAR_RUN_DIR="${SCRATCH:-$HOME}/hawkes_bivar_godambe_clean_$(date +%Y%m%d_%H%M%S)"
./run_quick_pilot_bivar.sh
./submit_slurm_writable_bivar.sh "$HAWKES_BIVAR_RUN_DIR"
```

After the Slurm combine job finishes:

```bash
./check_bivar_outputs.sh \
  "$HAWKES_BIVAR_RUN_DIR/bivar_results" \
  "$HAWKES_BIVAR_RUN_DIR/bivar_population" \
  "$HAWKES_BIVAR_RUN_DIR/bivar_figures" \
  "$HAWKES_BIVAR_RUN_DIR/bivar_figures_paper"
```

## Manuscript figure

Copy the side-by-side figure into your article's `plots/` folder:

```bash
mkdir -p /path/to/article/plots
cp "$HAWKES_BIVAR_RUN_DIR/bivar_figures_paper/paper_bivar_efficiency_pair_Tge1000.pdf" \
   /path/to/article/plots/
```

In the TeX file, include it as:

```latex
\begin{figure}[!t]
\centering
\includegraphics[width=\textwidth]{plots/paper_bivar_efficiency_pair_Tge1000}
\caption{Bivariate linear Hawkes simulation. Left: overall normalized RMSE, with the dotted line denoting MLE parity and horizontal lines denoting population Godambe targets. Right: population asymptotic standard-error inflation by parameter, relative to the MLE benchmark.}
\label{fig:bivar-overall-efficiency}
\label{fig:bivar-godambe-parameter}
\end{figure}
```

Or patch the old manuscript automatically:

```bash
python3 patch_tex_bivar_side_by_side.py 'main(3).tex' main_side_by_side.tex
```

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

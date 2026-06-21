# Terminal steps for a clean rerun

This package separates the code directory from the run directory. Keep the unpacked source directory clean, and write all generated CSVs, logs, and figures to a writable scratch/project directory.

## 1. Unpack and enter the package

```bash
unzip HawkesAnnalStudy_clean_final.zip
cd HawkesAnnalStudy_clean_final
chmod +x *.sh *.R *.py
```

Load R and make sure `Rcpp` is available:

```bash
module load R
Rscript -e 'if (!requireNamespace("Rcpp", quietly=TRUE)) install.packages("Rcpp", repos="https://cloud.r-project.org")'
```

## 2. Choose a writable run directory

```bash
export HAWKES_BIVAR_RUN_DIR="${SCRATCH:-$HOME}/hawkes_bivar_godambe_clean_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$HAWKES_BIVAR_RUN_DIR"
```

## 3. Run a local smoke test first

```bash
./run_quick_pilot_bivar.sh
```

The pilot also tests the final side-by-side article figure generator:

```text
$HAWKES_BIVAR_RUN_DIR/pilot_figures_paper/paper_bivar_pilot_efficiency_pair_Tge300.pdf
$HAWKES_BIVAR_RUN_DIR/pilot_figures_paper/paper_bivar_pilot_efficiency_pair_Tge300.png
```

## 4. Submit the full Slurm rerun

```bash
./submit_slurm_writable_bivar.sh "$HAWKES_BIVAR_RUN_DIR"
```

The wrapper submits the population array, estimator array, and dependent combine/figure job.

## 5. Check outputs after the combine job finishes

```bash
./check_bivar_outputs.sh \
  "$HAWKES_BIVAR_RUN_DIR/bivar_results" \
  "$HAWKES_BIVAR_RUN_DIR/bivar_population" \
  "$HAWKES_BIVAR_RUN_DIR/bivar_figures" \
  "$HAWKES_BIVAR_RUN_DIR/bivar_figures_paper"
```

The manuscript-ready side-by-side figure should be:

```bash
ls -lh "$HAWKES_BIVAR_RUN_DIR/bivar_figures_paper/paper_bivar_efficiency_pair_Tge1000.pdf"
```

## 6. Copy the figure into the article

```bash
mkdir -p /path/to/article/plots
cp "$HAWKES_BIVAR_RUN_DIR/bivar_figures_paper/paper_bivar_efficiency_pair_Tge1000.pdf" \
   /path/to/article/plots/
```

Use this in LaTeX:

```latex
\includegraphics[width=\textwidth]{plots/paper_bivar_efficiency_pair_Tge1000}
```

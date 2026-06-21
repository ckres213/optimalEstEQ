# Clean Slurm commands

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

After the dependent combine job finishes:

```bash
./check_bivar_outputs.sh \
  "$HAWKES_BIVAR_RUN_DIR/bivar_results" \
  "$HAWKES_BIVAR_RUN_DIR/bivar_population" \
  "$HAWKES_BIVAR_RUN_DIR/bivar_figures" \
  "$HAWKES_BIVAR_RUN_DIR/bivar_figures_paper"

ls -lh "$HAWKES_BIVAR_RUN_DIR/bivar_figures_paper/paper_bivar_efficiency_pair_Tge1000.pdf"
```

Copy the manuscript figure into the article:

```bash
mkdir -p /path/to/article/plots
cp "$HAWKES_BIVAR_RUN_DIR/bivar_figures_paper/paper_bivar_efficiency_pair_Tge1000.pdf" \
   /path/to/article/plots/
```

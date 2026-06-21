# Implementation notes

This bundle keeps the cleaned bivariate Hawkes Godambe simulation code and adds the final manuscript-side figure generator, `make_article_bivar_pair.R`.

The augmented estimator is the overidentified non-score GMM estimator already present in the cleaned study. The final manuscript figure is generated after combining the population and estimator arrays, using:

```bash
Rscript make_article_bivar_pair.R \
  --indir "$HAWKES_BIVAR_RUN_DIR/bivar_results" \
  --popdir "$HAWKES_BIVAR_RUN_DIR/bivar_population" \
  --outdir "$HAWKES_BIVAR_RUN_DIR/bivar_figures_paper" \
  --prefix paper_bivar \
  --tmin 1000
```

The output files are:

```text
paper_bivar_efficiency_pair_Tge1000.pdf
paper_bivar_efficiency_pair_Tge1000.png
```

The Slurm combine script and the quick pilot script both call this generator automatically.

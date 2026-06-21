# Article figure commands

After the full rerun finishes, copy the final side-by-side PDF into your manuscript's `plots/` directory:

```bash
ARTICLE_DIR="/path/to/your/article"
mkdir -p "$ARTICLE_DIR/plots"
cp "$HAWKES_BIVAR_RUN_DIR/bivar_figures_paper/paper_bivar_efficiency_pair_Tge1000.pdf" \
   "$ARTICLE_DIR/plots/"
```

Use the figure environment in `ARTICLE_FIGURE_SNIPPET.tex`, or patch the manuscript automatically:

```bash
python3 patch_tex_bivar_side_by_side.py "$ARTICLE_DIR/main(3).tex" \
  "$ARTICLE_DIR/main_side_by_side.tex"
```

Then compile:

```bash
cd "$ARTICLE_DIR"
latexmk -pdf main_side_by_side.tex
```

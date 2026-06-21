#!/usr/bin/env Rscript

## Full-width, two-panel article figure for the bivariate Hawkes study.
## Writes both PDF and PNG versions. The default layout is intended to be
## included in LaTeX as:
##   \includegraphics[width=\textwidth]{plots/paper_bivar_efficiency_pair_Tge1000}
##
## Normal clean-package usage:
##   Rscript make_article_bivar_pair.R \
##     --indir bivar_results \
##     --popdir bivar_population \
##     --outdir bivar_figures_paper \
##     --prefix paper_bivar \
##     --tmin 1000
##
## Backward-compatible usage for an old completed run directory:
##   Rscript make_article_bivar_pair.R hawkes_bivar_godambe_run1 1000

parse_args <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  out <- list(
    indir = "bivar_results",
    popdir = NA_character_,
    outdir = "bivar_figures_paper",
    prefix = "paper_bivar",
    tmin = "1000",
    pdf_width = "7.6",
    pdf_height = "4.5",
    png_width = "2600",
    png_height = "1550",
    png_res = "300",
    pointsize = "14"
  )

  ## Positional mode lets the script be dropped into the original old study
  ## directory and run directly against hawkes_bivar_godambe_run1.
  if (length(a) >= 1L && !grepl("^--", a[[1L]])) {
    run_dir <- a[[1L]]
    out$indir <- file.path(run_dir, "bivar_results")
    out$popdir <- file.path(run_dir, "bivar_population")
    out$outdir <- file.path(run_dir, "bivar_figures_paper")
    if (length(a) >= 2L) out$tmin <- a[[2L]]
    if (length(a) > 2L) {
      stop("Too many positional arguments. Use: Rscript make_article_bivar_pair.R RUN_DIR [TMIN]", call. = FALSE)
    }
  } else {
    i <- 1L
    while (i <= length(a)) {
      if (!grepl("^--", a[[i]])) stop("Expected a --key argument, got: ", a[[i]], call. = FALSE)
      key <- gsub("-", "_", sub("^--", "", a[[i]]))
      if (i == length(a)) stop("Missing value for --", key, call. = FALSE)
      out[[key]] <- a[[i + 1L]]
      i <- i + 2L
    }
  }

  out$tmin <- as.numeric(out$tmin)
  out$pdf_width <- as.numeric(out$pdf_width)
  out$pdf_height <- as.numeric(out$pdf_height)
  out$png_width <- as.integer(out$png_width)
  out$png_height <- as.integer(out$png_height)
  out$png_res <- as.integer(out$png_res)
  out$pointsize <- as.numeric(out$pointsize)

  if (!is.finite(out$tmin)) stop("tmin must be numeric", call. = FALSE)
  if (!is.finite(out$pdf_width) || !is.finite(out$pdf_height)) stop("PDF dimensions must be numeric", call. = FALSE)
  if (!is.finite(out$pointsize)) stop("pointsize must be numeric", call. = FALSE)
  out
}

stopf <- function(...) stop(paste0(...), call. = FALSE)

read_csv_required <- function(path) {
  if (!file.exists(path)) stopf("Missing file: ", path)
  if (file.info(path)$size <= 0) stopf("Empty file: ", path)
  read.csv(path, stringsAsFactors = FALSE)
}

find_input <- function(filename, indir, popdir = NA_character_, required = TRUE) {
  candidates <- c(file.path(indir, filename))
  if (!is.na(popdir) && nzchar(popdir)) {
    candidates <- c(candidates, file.path(popdir, filename))
  }
  hit <- candidates[file.exists(candidates) & file.info(candidates)$size > 0]
  if (length(hit)) return(hit[1L])
  if (required) stopf("Could not find ", filename, " in indir or popdir")
  NA_character_
}

safe_range <- function(x, include = NULL, pad = 0.08, lower_floor = NULL) {
  x <- c(x, include)
  x <- x[is.finite(x)]
  if (!length(x)) return(c(0, 1))
  r <- range(x)
  if (diff(r) == 0) {
    dd <- max(0.05, 0.1 * abs(r[1]))
    r <- r + c(-dd, dd)
  }
  out <- r + c(-pad, pad) * diff(r)
  if (!is.null(lower_floor)) out[1L] <- max(lower_floor, out[1L])
  out
}

ensure_group <- function(d) {
  if ("group" %in% names(d)) return(d)
  if (!"param" %in% names(d)) return(d)
  g <- rep("other", nrow(d))
  g[d$param %in% c("mu1", "mu2")] <- "baselines"
  g[d$param %in% c("alpha11", "alpha22")] <- "self-excitation"
  g[d$param %in% c("alpha12", "alpha21")] <- "cross-excitation"
  g[d$param %in% c("beta")] <- "decay"
  d$group <- g
  d
}

param_order <- c("mu1", "mu2", "alpha11", "alpha22", "alpha12", "alpha21", "beta")

param_plot_labels <- function(params) {
  labs <- lapply(params, function(x) {
    switch(
      x,
      mu1 = expression(mu[1])[[1]],
      mu2 = expression(mu[2])[[1]],
      alpha11 = expression(alpha[11])[[1]],
      alpha22 = expression(alpha[22])[[1]],
      alpha12 = expression(alpha[12])[[1]],
      alpha21 = expression(alpha[21])[[1]],
      beta = expression(beta)[[1]],
      as.name(x)
    )
  })
  as.expression(labs)
}

method_cols <- c(
  MLE = "#1f77b4",
  GMM_Dtheta = "#d95f02",
  GMM_aug = "#1b9e77"
)

method_pch <- c(
  MLE = 16,
  GMM_Dtheta = 17,
  GMM_aug = 15
)

method_lty <- c(
  MLE = 1,
  GMM_Dtheta = 2,
  GMM_aug = 4
)

pretty_method <- function(m) {
  out <- m
  out[m == "MLE"] <- "MLE"
  out[m == "GMM_Dtheta"] <- "just-id GMM"
  out[m == "GMM_aug"] <- "over-id GMM"
  out
}

method_order <- function(x) {
  match(x, c("MLE", "GMM_Dtheta", "GMM_aug"))
}

aggregate_overall_efficiency <- function(indir, popdir) {
  std_path <- find_input("bivar_standardized_summary.csv", indir, popdir, required = TRUE)
  pop_path <- find_input("population_parameter_inflation.csv", indir, popdir, required = TRUE)

  std <- read_csv_required(std_path)
  pop <- read_csv_required(pop_path)

  std <- ensure_group(std)
  pop <- ensure_group(pop)

  need_std <- c("T", "method", "param", "scaled_rmse_over_target")
  need_pop <- c("param", "se_inflation_gmm_to_mle")
  miss_std <- setdiff(need_std, names(std))
  miss_pop <- setdiff(need_pop, names(pop))
  if (length(miss_std)) stopf("Missing columns in bivar_standardized_summary.csv: ", paste(miss_std, collapse = ", "))
  if (length(miss_pop)) stopf("Missing columns in population_parameter_inflation.csv: ", paste(miss_pop, collapse = ", "))

  std$T <- as.numeric(std$T)
  std$scaled_rmse_over_target <- as.numeric(std$scaled_rmse_over_target)
  pop$se_inflation_gmm_to_mle <- as.numeric(pop$se_inflation_gmm_to_mle)
  if ("se_inflation_aug_to_mle" %in% names(pop)) {
    pop$se_inflation_aug_to_mle <- as.numeric(pop$se_inflation_aug_to_mle)
  }

  keep_pop <- intersect(
    c("param", "se_inflation_gmm_to_mle", "se_inflation_aug_to_mle"),
    names(pop)
  )

  d <- merge(std, pop[, keep_pop, drop = FALSE], by = "param", all.x = TRUE)

  d$se_inflation_to_mle <- 1
  idx <- d$method == "GMM_Dtheta"
  d$se_inflation_to_mle[idx] <- d$se_inflation_gmm_to_mle[idx]

  if ("se_inflation_aug_to_mle" %in% names(d)) {
    idx <- d$method == "GMM_aug"
    d$se_inflation_to_mle[idx] <- d$se_inflation_aug_to_mle[idx]
  }

  d$efficiency_rmse_mle_scale <- d$scaled_rmse_over_target * d$se_inflation_to_mle

  split_key <- paste(d$T, d$method, sep = "|")
  overall <- do.call(rbind, lapply(split(d, split_key), function(z) {
    data.frame(
      T = z$T[1L],
      method = z$method[1L],
      n_param = sum(is.finite(z$efficiency_rmse_mle_scale)),
      overall_efficiency_rmse_mle_scale =
        sqrt(mean(z$efficiency_rmse_mle_scale^2, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }))
  rownames(overall) <- NULL
  overall <- overall[order(overall$T, method_order(overall$method)), , drop = FALSE]

  pop_target <- data.frame(
    method = "MLE",
    population_target = 1,
    stringsAsFactors = FALSE
  )

  pop_target <- rbind(
    pop_target,
    data.frame(
      method = "GMM_Dtheta",
      population_target = sqrt(mean(pop$se_inflation_gmm_to_mle^2, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  )

  if ("se_inflation_aug_to_mle" %in% names(pop) &&
      any(is.finite(pop$se_inflation_aug_to_mle))) {
    pop_target <- rbind(
      pop_target,
      data.frame(
        method = "GMM_aug",
        population_target = sqrt(mean(pop$se_inflation_aug_to_mle^2, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
    )
  }

  list(overall = overall, pop_target = pop_target, pop = pop)
}

make_parameter_table <- function(indir, popdir) {
  pop_path <- find_input("population_parameter_inflation.csv", indir, popdir, required = TRUE)
  parinf <- read_csv_required(pop_path)

  need <- c("param", "se_inflation_gmm_to_mle")
  miss <- setdiff(need, names(parinf))
  if (length(miss)) stopf("Missing columns in population_parameter_inflation.csv: ", paste(miss, collapse = ", "))

  parinf$param <- factor(parinf$param, levels = param_order)
  parinf <- parinf[order(parinf$param), , drop = FALSE]
  parinf <- parinf[!is.na(parinf$param), , drop = FALSE]

  parinf$se_inflation_gmm_to_mle <- as.numeric(parinf$se_inflation_gmm_to_mle)
  if ("se_inflation_aug_to_mle" %in% names(parinf)) {
    parinf$se_inflation_aug_to_mle <- as.numeric(parinf$se_inflation_aug_to_mle)
  }

  par_long <- data.frame(
    method = "GMM_Dtheta",
    param = as.character(parinf$param),
    se_ratio = parinf$se_inflation_gmm_to_mle,
    stringsAsFactors = FALSE
  )

  if ("se_inflation_aug_to_mle" %in% names(parinf) &&
      any(is.finite(parinf$se_inflation_aug_to_mle))) {
    par_long <- rbind(
      par_long,
      data.frame(
        method = "GMM_aug",
        param = as.character(parinf$param),
        se_ratio = parinf$se_inflation_aug_to_mle,
        stringsAsFactors = FALSE
      )
    )
  }

  par_long <- par_long[is.finite(par_long$se_ratio), , drop = FALSE]
  list(parinf = parinf, par_long = par_long)
}

draw_article_pair <- function(indir, popdir, tmin) {
  tabs <- aggregate_overall_efficiency(indir, popdir)
  d <- tabs$overall[tabs$overall$T >= tmin, , drop = FALSE]
  pop_target <- tabs$pop_target
  if (!nrow(d)) stopf("No overall rows with T >= ", tmin)

  ptabs <- make_parameter_table(indir, popdir)
  parinf <- ptabs$parinf
  par_long <- ptabs$par_long
  gmm_methods <- c("GMM_Dtheta", "GMM_aug")
  gmm_methods <- gmm_methods[gmm_methods %in% unique(par_long$method)]

  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)

  ## Top row: the two scientific panels. Bottom row: one shared legend strip.
  layout(
    matrix(c(1, 2,
             3, 3), nrow = 2, byrow = TRUE),
    widths = c(1.18, 0.98),
    heights = c(1.00, 0.24)
  )
  par(oma = c(0, 0, 0, 0), mgp = c(2.55, 0.78, 0), tcl = -0.25)

  ## Panel A: overall efficiency/RMSE comparison.
  par(mar = c(4.7, 5.0, 1.9, 1.1))
  yy <- c(d$overall_efficiency_rmse_mle_scale, pop_target$population_target, 1)
  plot(
    NA,
    xlim = range(d$T),
    ylim = safe_range(yy, include = 1, pad = 0.12),
    log = "x",
    xaxt = "n",
    xlab = "T",
    ylab = "overall normalized RMSE",
    main = "",
    cex.lab = 1.22,
    cex.axis = 1.04
  )
  ticks <- sort(unique(d$T))
  axis(1, at = ticks, labels = format(ticks, scientific = FALSE), cex.axis = 0.98)
  grid(col = "grey90")
  abline(h = 1, lty = 3, lwd = 1.4, col = "grey40")

  methods <- names(method_cols)[names(method_cols) %in% unique(d$method)]
  for (m in methods) {
    z <- d[d$method == m, , drop = FALSE]
    z <- z[order(z$T), , drop = FALSE]
    lines(
      z$T,
      z$overall_efficiency_rmse_mle_scale,
      type = "b",
      pch = method_pch[[m]],
      lty = method_lty[[m]],
      col = method_cols[[m]],
      lwd = 2.4,
      cex = 1.05
    )
    target <- pop_target$population_target[pop_target$method == m]
    if (length(target) && is.finite(target[1L])) {
      abline(
        h = target[1L],
        col = adjustcolor(method_cols[[m]], alpha.f = 0.40),
        lwd = 1.8,
        lty = 1
      )
    }
  }
  mtext("(a)", side = 3, adj = 0, line = 0.35, font = 2, cex = 1.18)

  ## Panel B: coordinatewise population Godambe inflation.
  par(mar = c(4.7, 6.3, 1.9, 0.8))
  params <- as.character(parinf$param)
  y <- seq_along(params)
  offsets <- if (length(gmm_methods) == 1L) {
    setNames(0, gmm_methods)
  } else {
    setNames(seq(-0.14, 0.14, length.out = length(gmm_methods)), gmm_methods)
  }
  xr <- safe_range(par_long$se_ratio, include = 1, pad = 0.12)
  xr[1L] <- min(0.98, xr[1L])
  plot(
    NA,
    xlim = xr,
    ylim = c(0.4, length(params) + 0.8),
    yaxt = "n",
    xlab = "asymptotic SE ratio / MLE",
    ylab = "",
    main = "",
    cex.lab = 1.22,
    cex.axis = 1.04
  )
  axis(2, at = y, labels = param_plot_labels(params), las = 2, cex.axis = 1.08)
  grid(nx = NULL, ny = NA, col = "grey90")
  abline(v = 1, lty = 3, lwd = 1.4, col = "grey40")

  for (m in gmm_methods) {
    z <- par_long[par_long$method == m, , drop = FALSE]
    yy <- y[match(z$param, params)] + offsets[[m]]
    segments(
      1,
      yy,
      z$se_ratio,
      yy,
      col = method_cols[[m]],
      lwd = 4.0,
      lend = "butt"
    )
    points(
      z$se_ratio,
      yy,
      pch = method_pch[[m]],
      col = method_cols[[m]],
      cex = 0.95
    )
  }
  mtext("(b)", side = 3, adj = 0, line = 0.35, font = 2, cex = 1.18)

  ## Shared bottom legend strip. This keeps both plotting regions wide and
  ## avoids drawing text over the data.
  par(mar = c(0.2, 0.2, 0.1, 0.2))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))

  legend(
    x = 0.28,
    y = 0.58,
    legend = pretty_method(methods),
    col = method_cols[methods],
    pch = method_pch[methods],
    lty = method_lty[methods],
    lwd = 2.4,
    bty = "n",
    horiz = TRUE,
    xjust = 0.5,
    yjust = 0.5,
    cex = 0.98,
    pt.cex = 1.05,
    seg.len = 2.8,
    xpd = NA
  )

  legend(
    x = 0.78,
    y = 0.58,
    legend = c("population target", "MLE parity"),
    col = c("grey45", "grey45"),
    lty = c(1, 3),
    lwd = c(1.8, 1.4),
    bty = "n",
    horiz = TRUE,
    xjust = 0.5,
    yjust = 0.5,
    cex = 0.98,
    seg.len = 2.8,
    xpd = NA
  )
}

main <- function() {
  args <- parse_args()
  dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

  stem <- file.path(args$outdir, paste0(args$prefix, "_efficiency_pair_Tge", args$tmin))
  out_pdf <- paste0(stem, ".pdf")
  out_png <- paste0(stem, ".png")

  pdf(
    out_pdf,
    width = args$pdf_width,
    height = args$pdf_height,
    onefile = FALSE,
    pointsize = args$pointsize,
    useDingbats = FALSE
  )
  draw_article_pair(args$indir, args$popdir, args$tmin)
  dev.off()

  png(
    out_png,
    width = args$png_width,
    height = args$png_height,
    res = args$png_res,
    pointsize = args$pointsize
  )
  draw_article_pair(args$indir, args$popdir, args$tmin)
  dev.off()

  cat("Wrote article side-by-side figure:\n")
  cat("  ", out_pdf, "\n", sep = "")
  cat("  ", out_png, "\n", sep = "")
}

main()

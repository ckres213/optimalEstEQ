#!/usr/bin/env Rscript

parse_args <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  out <- list(
    indir = "bivar_results",
    popdir = NA_character_,
    outdir = "bivar_figures_paper",
    prefix = "paper_bivar",
    tmin = "1000"
  )
  i <- 1L
  while (i <= length(a)) {
    key <- sub("^--", "", a[[i]])
    if (i == length(a)) stop("Missing value for --", key, call. = FALSE)
    out[[key]] <- a[[i + 1L]]
    i <- i + 2L
  }
  out$tmin <- as.numeric(out$tmin)
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
  if (length(hit)) return(hit[1])
  if (required) stopf("Could not find ", filename, " in indir or popdir")
  NA_character_
}

save_png <- function(path, width = 1700, height = 1050, res = 170, draw_fun) {
  png(path, width = width, height = height, res = res)
  draw_fun()
  dev.off()
}

save_pdf <- function(path, width = 7.2, height = 5.8, draw_fun) {
  pdf(path, width = width, height = height, onefile = FALSE)
  draw_fun()
  dev.off()
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
  if (!is.null(lower_floor)) out[1] <- max(lower_floor, out[1])
  out
}

log_x_axis <- function(T) {
  ticks <- sort(unique(as.numeric(T)))
  axis(1, at = ticks, labels = ticks)
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
  MLE = "#0072B2",
  GMM_Dtheta = "#D55E00",
  GMM_aug = "#009E73"
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
      T = z$T[1],
      method = z$method[1],
      n_param = sum(is.finite(z$efficiency_rmse_mle_scale)),
      overall_efficiency_rmse_mle_scale =
        sqrt(mean(z$efficiency_rmse_mle_scale^2, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }))
  rownames(overall) <- NULL
  overall <- overall[order(overall$T, method_order(overall$method)), ]

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

  list(overall = overall, pop_target = pop_target)
}

plot_overall_efficiency <- function(indir, popdir, outdir, prefix, tmin) {
  tabs <- aggregate_overall_efficiency(indir, popdir)
  d <- tabs$overall[tabs$overall$T >= tmin, , drop = FALSE]
  pop_target <- tabs$pop_target
  if (!nrow(d)) stopf("No overall rows with T >= ", tmin)

  out_png <- file.path(outdir, paste0(prefix, "_overall_efficiency_Tge", tmin, ".png"))

  draw_fun <- function() {
    old <- par(no.readonly = TRUE)
    on.exit(par(old), add = TRUE)

    layout(
      matrix(c(1, 2), nrow = 2),
      heights = c(1.0, 0.16)
    )

    par(mar = c(5.0, 5.5, 1.0, 1.2))

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
      cex.lab = 1.15,
      cex.axis = 1.05
    )

    log_x_axis(d$T)
    grid(col = "grey88")
    abline(h = 1, lty = 3, lwd = 1.5, col = "grey45")

    methods <- intersect(names(method_cols), unique(d$method))
    for (m in methods) {
      z <- d[d$method == m, , drop = FALSE]
      z <- z[order(z$T), ]

      lines(
        z$T,
        z$overall_efficiency_rmse_mle_scale,
        type = "b",
        pch = method_pch[[m]],
        lty = method_lty[[m]],
        col = method_cols[[m]],
        lwd = 3.2,
        cex = 1.25
      )

      target <- pop_target$population_target[pop_target$method == m]
      if (length(target) && is.finite(target[1])) {
        abline(h = target[1], col = adjustcolor(method_cols[[m]], alpha.f = 0.45), lwd = 3)
      }
    }

    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend(
      "center",
      legend = c(paste0(pretty_method(methods), " empirical"),
                 "colored horizontal lines = population targets",
                 "dotted horizontal line = parity"),
      col = c(method_cols[methods], "grey45", "grey25"),
      pch = c(method_pch[methods], NA, NA),
      lty = c(method_lty[methods], 1, 3),
      lwd = c(rep(3.2, length(methods)), 3, 1.7),
      bty = "n",
      cex = 0.88,
      ncol = 2
    )
  }

  save_png(out_png, width = 1700, height = 1050, res = 170, draw_fun = draw_fun)
  save_pdf(sub("\\.png$", ".pdf", out_png), width = 8.6, height = 5.4, draw_fun = draw_fun)

  out_png
}

plot_per_parameter_godambe <- function(indir, popdir, outdir, prefix) {
  pop_path <- find_input("population_parameter_inflation.csv", indir, popdir, required = TRUE)
  parinf <- read_csv_required(pop_path)

  need <- c("param", "se_inflation_gmm_to_mle")
  miss <- setdiff(need, names(parinf))
  if (length(miss)) stopf("Missing columns in population_parameter_inflation.csv: ", paste(miss, collapse = ", "))

  parinf$param <- factor(parinf$param, levels = param_order)
  parinf <- parinf[order(parinf$param), ]
  parinf <- parinf[!is.na(parinf$param), ]

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
  gmm_methods <- intersect(c("GMM_Dtheta", "GMM_aug"), unique(par_long$method))

  out_png <- file.path(outdir, paste0(prefix, "_per_parameter_godambe_se_inflation.png"))

  draw_fun <- function() {
    old <- par(no.readonly = TRUE)
    on.exit(par(old), add = TRUE)

    layout(
      matrix(c(1, 2), nrow = 2),
      heights = c(1.0, 0.15)
    )

    par(mar = c(5.1, 7.8, 1.0, 1.2))

    params <- as.character(parinf$param)
    y <- seq_along(params)

    offsets <- if (length(gmm_methods) == 1L) {
      setNames(0, gmm_methods)
    } else {
      setNames(seq(-0.15, 0.15, length.out = length(gmm_methods)), gmm_methods)
    }

    xr <- safe_range(par_long$se_ratio, include = 1, pad = 0.14)
    xr[1] <- min(0.98, xr[1])

    plot(
      NA,
      xlim = xr,
      ylim = c(0.4, length(params) + 0.6),
      yaxt = "n",
      xlab = "asymptotic SE ratio relative to MLE",
      ylab = "",
      main = "",
      cex.lab = 1.12,
      cex.axis = 1.02
    )

    axis(2, at = y, labels = param_plot_labels(params), las = 2, cex.axis = 1.02)
    grid(nx = NULL, ny = NA, col = "grey88")
    abline(v = 1, lty = 3, lwd = 1.5, col = "grey45")

    for (m in gmm_methods) {
      z <- par_long[par_long$method == m, , drop = FALSE]
      yy <- y[match(z$param, params)] + offsets[[m]]
      segments(1, yy, z$se_ratio, yy, col = method_cols[[m]], lwd = 5, lend = "butt")
      # Endpoint symbols intentionally omitted; horizontal segments carry the comparison.
    }

    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend(
      "center",
      legend = c(pretty_method(gmm_methods), "MLE benchmark"),
      col = c(method_cols[gmm_methods], "grey25"),
      pch = rep(NA, length(gmm_methods) + 1),
      lty = c(method_lty[gmm_methods], 3),
      lwd = c(rep(5, length(gmm_methods)), 1.8),
      bty = "n",
      cex = 0.84,
      ncol = length(gmm_methods) + 1
    )
  }

  save_png(out_png, width = 1500, height = 1050, res = 170, draw_fun = draw_fun)
  save_pdf(sub("\\.png$", ".pdf", out_png), width = 7.2, height = 5.4, draw_fun = draw_fun)

  out_png
}

load_ci_overall <- function(indir, popdir) {
  byp_path <- find_input("bivar_asymptotic_ci_by_parameter.csv", indir, popdir, required = FALSE)
  overall_path <- find_input("bivar_asymptotic_ci_overall.csv", indir, popdir, required = FALSE)

  if (!is.na(byp_path)) {
    d <- read_csv_required(byp_path)
    need <- c("T", "method", "coverage", "ci_level")
    miss <- setdiff(need, names(d))
    if (length(miss)) stopf("Missing columns in bivar_asymptotic_ci_by_parameter.csv: ", paste(miss, collapse = ", "))

    d$T <- as.numeric(d$T)
    d$coverage <- as.numeric(d$coverage)
    d$ci_level <- as.numeric(d$ci_level)

    if ("mean_ci_width" %in% names(d)) {
      d$ci_width_for_plot <- as.numeric(d$mean_ci_width)
    } else if ("ci_width" %in% names(d)) {
      d$ci_width_for_plot <- as.numeric(d$ci_width)
    } else if ("ci_half_width" %in% names(d)) {
      d$ci_width_for_plot <- 2 * as.numeric(d$ci_half_width)
    } else if ("asymptotic_se" %in% names(d) && "zcrit" %in% names(d)) {
      d$ci_width_for_plot <- 2 * as.numeric(d$zcrit) * as.numeric(d$asymptotic_se)
    } else {
      stopf("Could not find a CI width column in bivar_asymptotic_ci_by_parameter.csv. Available columns: ",
            paste(names(d), collapse = ", "))
    }

    key <- paste(d$T, d$method, sep = "|")
    out <- do.call(rbind, lapply(split(d, key), function(z) {
      data.frame(
        T = z$T[1],
        method = z$method[1],
        ci_level = z$ci_level[is.finite(z$ci_level)][1],
        coverage = mean(z$coverage, na.rm = TRUE),
        mean_ci_width = mean(z$ci_width_for_plot, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
    rownames(out) <- NULL
    out <- out[order(out$T, method_order(out$method)), ]
    return(out)
  }

  if (!is.na(overall_path)) {
    d <- read_csv_required(overall_path)
    need <- c("T", "method", "coverage", "ci_level")
    miss <- setdiff(need, names(d))
    if (length(miss)) stopf("Missing columns in bivar_asymptotic_ci_overall.csv: ", paste(miss, collapse = ", "))
    d$T <- as.numeric(d$T)
    d$coverage <- as.numeric(d$coverage)
    d$ci_level <- as.numeric(d$ci_level)
    if ("mean_ci_width" %in% names(d)) {
      d$mean_ci_width <- as.numeric(d$mean_ci_width)
    } else {
      stopf("Overall CI file lacks mean_ci_width, and by-parameter file was not found.")
    }
    return(d)
  }

  stopf("Could not find CI summary files.")
}

plot_appendix_ci <- function(indir, popdir, outdir, prefix, tmin) {
  d <- load_ci_overall(indir, popdir)
  d <- d[is.finite(d$T) & d$T >= tmin, , drop = FALSE]
  if (!nrow(d)) stopf("No CI rows with T >= ", tmin)

  target <- d$ci_level[is.finite(d$ci_level)][1]
  if (!is.finite(target)) target <- 0.95

  out_png <- file.path(outdir, paste0(prefix, "_appendix_ci_coverage_width_Tge", tmin, ".png"))

  draw_fun <- function() {
    old <- par(no.readonly = TRUE)
    on.exit(par(old), add = TRUE)

    layout(
      matrix(c(1, 2, 3), nrow = 1),
      widths = c(1.0, 1.0, 0.26)
    )

    methods <- intersect(names(method_cols), unique(d$method))

    # --------------------------
    # Left panel: CI coverage
    # --------------------------
    par(mar = c(5.0, 5.3, 1.0, 1.1))

    plot(
      NA,
      xlim = range(d$T),
      ylim = c(0.90, 1.00),
      log = "x",
      xaxt = "n",
      xlab = "T",
      ylab = "empirical coverage",
      main = "",
      cex.lab = 1.08,
      cex.axis = 0.98
    )

    log_x_axis(d$T)
    grid(col = "grey88")
    abline(h = target, lty = 3, lwd = 1.5, col = "grey45")

    for (m in methods) {
      z <- d[d$method == m, , drop = FALSE]
      z <- z[order(z$T), ]
      lines(
        z$T,
        z$coverage,
        type = "b",
        pch = method_pch[[m]],
        lty = method_lty[[m]],
        col = method_cols[[m]],
        lwd = 2.8,
        cex = 1.1
      )
    }

    # --------------------------
    # Middle panel: CI width
    # --------------------------
    par(mar = c(5.0, 5.3, 1.0, 1.1))

    yy <- d$mean_ci_width

    plot(
      NA,
      xlim = range(d$T),
      ylim = safe_range(yy, pad = 0.08),
      log = "x",
      xaxt = "n",
      xlab = "T",
      ylab = "mean CI width",
      main = "",
      cex.lab = 1.08,
      cex.axis = 0.98
    )

    log_x_axis(d$T)
    grid(col = "grey88")

    for (m in methods) {
      z <- d[d$method == m, , drop = FALSE]
      z <- z[order(z$T), ]
      lines(
        z$T,
        z$mean_ci_width,
        type = "b",
        pch = method_pch[[m]],
        lty = method_lty[[m]],
        col = method_cols[[m]],
        lwd = 2.8,
        cex = 1.1
      )
    }

    # --------------------------
    # Right column: vertical legend
    # --------------------------
    par(mar = c(0, 0, 0, 0))
    plot.new()

    legend(
      "center",
      legend = c(pretty_method(methods), paste0("target = ", sprintf("%.2f", target))),
      col = c(method_cols[methods], "grey45"),
      pch = c(method_pch[methods], NA),
      lty = c(method_lty[methods], 3),
      lwd = c(rep(2.8, length(methods)), 1.5),
      bty = "n",
      cex = 0.84,
      ncol = 1,
      y.intersp = 1.35,
      x.intersp = 0.75,
      seg.len = 1.5
    )
  }

  save_png(out_png, width = 2150, height = 1050, res = 170, draw_fun = draw_fun)
  save_pdf(sub("\\.png$", ".pdf", out_png), width = 10.6, height = 5.2, draw_fun = draw_fun)

  out_png
}


plot_covariance_eigenvalues <- function(indir, popdir, outdir, prefix) {
  eig_path <- find_input("population_relative_eigenvalues.csv", indir, popdir, required = TRUE)
  eig <- read_csv_required(eig_path)

  need <- c("eig_index", "covariance_eigenvalue")
  miss <- setdiff(need, names(eig))
  if (length(miss)) stopf("Missing columns in population_relative_eigenvalues.csv: ", paste(miss, collapse = ", "))

  eig$eig_index <- as.numeric(eig$eig_index)
  eig$covariance_eigenvalue <- as.numeric(eig$covariance_eigenvalue)
  if (!"method" %in% names(eig)) eig$method <- "GMM_Dtheta"

  eig <- eig[is.finite(eig$eig_index) & is.finite(eig$covariance_eigenvalue), , drop = FALSE]
  eig_methods <- intersect(c("GMM_Dtheta", "GMM_aug"), unique(eig$method))

  out_png <- file.path(outdir, paste0(prefix, "_covariance_eigenvalue_inflation.png"))

  draw_fun <- function() {
    old <- par(no.readonly = TRUE)
    on.exit(par(old), add = TRUE)

    layout(
      matrix(c(1, 2), nrow = 2),
      heights = c(1.0, 0.15)
    )

    par(mar = c(5.0, 5.5, 1.0, 1.2))
    yr <- safe_range(eig$covariance_eigenvalue, include = 1, pad = 0.12)

    plot(
      NA,
      xlim = range(eig$eig_index),
      ylim = yr,
      xlab = "ordered eigenvalue index",
      ylab = "eigenvalue ratio relative to MLE",
      main = "",
      cex.lab = 1.12,
      cex.axis = 1.02
    )

    grid(col = "grey88")
    abline(h = 1, lty = 3, lwd = 1.5, col = "grey45")

    for (m in eig_methods) {
      z <- eig[eig$method == m, , drop = FALSE]
      z <- z[order(z$eig_index), ]
      lines(
        z$eig_index,
        z$covariance_eigenvalue,
        type = "b",
        pch = method_pch[[m]],
        lty = method_lty[[m]],
        lwd = 3.0,
        col = method_cols[[m]]
      )
    }

    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend(
      "center",
      legend = c(pretty_method(eig_methods), "MLE benchmark"),
      col = c(method_cols[eig_methods], "grey25"),
      pch = c(method_pch[eig_methods], NA),
      lty = c(method_lty[eig_methods], 3),
      lwd = c(rep(3.0, length(eig_methods)), 1.8),
      bty = "n",
      cex = 0.84,
      ncol = length(eig_methods) + 1
    )
  }

  save_png(out_png, width = 1450, height = 1000, res = 170, draw_fun = draw_fun)
  save_pdf(sub("\\.png$", ".pdf", out_png), width = 7.0, height = 5.0, draw_fun = draw_fun)

  out_png
}

args <- parse_args()
dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)

outputs <- c(
  plot_overall_efficiency(args$indir, args$popdir, args$outdir, args$prefix, args$tmin),
  plot_per_parameter_godambe(args$indir, args$popdir, args$outdir, args$prefix),
  plot_appendix_ci(args$indir, args$popdir, args$outdir, args$prefix, args$tmin),
  plot_covariance_eigenvalues(args$indir, args$popdir, args$outdir, args$prefix)
)

cat("Wrote final figures:\n")
cat(paste0("  ", outputs, collapse = "\n"), "\n")

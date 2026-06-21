#!/usr/bin/env Rscript
################################################################################
# paper_bivar_four_plots.R
#
# Publication-style plots for the bivariate Hawkes Godambe study.
# Outputs PNGs only:
#   1. <prefix>_population_godambe_clear.png
#   2. <prefix>_overall_efficiency_Tge<tmin>.png
#   3. <prefix>_overall_calibration_Tge<tmin>.png
#   4. <prefix>_ci_coverage_width_Tge<tmin>.png
#   5. <prefix>_loglog_rmse.png
#
# Requires only base R.
################################################################################

parse_args <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  out <- list(
    indir = "bivar_results",
    popdir = NA_character_,
    outdir = "bivar_figures_paper",
    prefix = "paper",
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

read_csv_optional <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0) return(NULL)
  tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
}

find_input <- function(filename, indir, popdir = NA_character_, required = TRUE) {
  candidates <- c(file.path(indir, filename))
  if (!is.na(popdir) && nzchar(popdir)) candidates <- c(candidates, file.path(popdir, filename))
  hit <- candidates[file.exists(candidates) & file.info(candidates)$size > 0]
  if (length(hit)) return(hit[1])
  if (required) stopf("Could not find ", filename, " in indir or popdir")
  NA_character_
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
group_order <- c("baselines", "self-excitation", "cross-excitation", "decay")

method_cols <- c(MLE = "#0072B2", GMM_Dtheta = "#D55E00", GMM_aug = "#009E73")
method_pch <- c(MLE = 16, GMM_Dtheta = 17, GMM_aug = 15)
method_lty <- c(MLE = 1, GMM_Dtheta = 2, GMM_aug = 4)
param_cols <- c(
  mu1 = "#1b9e77",
  mu2 = "#d95f02",
  alpha11 = "#7570b3",
  alpha22 = "#e7298a",
  alpha12 = "#66a61e",
  alpha21 = "#e6ab02",
  beta = "#1f78b4"
)
param_pch <- c(mu1 = 16, mu2 = 17, alpha11 = 15, alpha22 = 18, alpha12 = 1, alpha21 = 2, beta = 0)
param_lty <- c(mu1 = 1, mu2 = 2, alpha11 = 1, alpha22 = 2, alpha12 = 1, alpha21 = 2, beta = 1)

pretty_method <- function(m) {
  out <- m
  out[m == "GMM_Dtheta"] <- "just-identified Dtheta GMM"
  out[m == "GMM_aug"] <- "overidentified bounded-inverse GMM"
  out
}

safe_range <- function(x, include = NULL, pad = 0.08, lower_floor = NULL) {
  x <- c(x, include)
  x <- x[is.finite(x)]
  if (!length(x)) return(c(0, 1))
  r <- range(x)
  if (!is.null(lower_floor)) r[1] <- max(lower_floor, r[1])
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

save_png <- function(path, width = 1800, height = 1150, res = 170, draw_fun) {
  png(path, width = width, height = height, res = res)
  draw_fun()
  dev.off()
}

make_overall_tables <- function(indir, popdir) {
  std_path <- find_input("bivar_standardized_summary.csv", indir, popdir, required = TRUE)
  pop_path <- find_input("population_parameter_inflation.csv", indir, popdir, required = TRUE)

  std <- read_csv_required(std_path)
  pop <- read_csv_required(pop_path)

  std <- ensure_group(std)
  pop <- ensure_group(pop)

  need_std <- c("T", "method", "param", "group", "scaled_rmse_over_target")
  need_pop <- c("param", "se_inflation_gmm_to_mle")

  miss_std <- setdiff(need_std, names(std))
  miss_pop <- setdiff(need_pop, names(pop))
  if (length(miss_std)) stopf("Missing columns in bivar_standardized_summary.csv: ", paste(miss_std, collapse = ", "))
  if (length(miss_pop)) stopf("Missing columns in population_parameter_inflation.csv: ", paste(miss_pop, collapse = ", "))

  std$T <- as.numeric(std$T)
  std$scaled_rmse_over_target <- as.numeric(std$scaled_rmse_over_target)
  pop$se_inflation_gmm_to_mle <- as.numeric(pop$se_inflation_gmm_to_mle)
  if ("se_inflation_aug_to_mle" %in% names(pop)) pop$se_inflation_aug_to_mle <- as.numeric(pop$se_inflation_aug_to_mle)

  keep_pop <- intersect(c("param", "se_inflation_gmm_to_mle", "se_inflation_aug_to_mle"), names(pop))
  d <- merge(std, pop[, keep_pop, drop = FALSE], by = "param", all.x = TRUE)
  d <- ensure_group(d)

  # Calibration scale:
  #   sqrt(T) * RMSE(method) / asymptotic SD(method).
  d$calibration_rmse <- d$scaled_rmse_over_target

  # Common MLE scale:
  #   sqrt(T) * RMSE(method) / asymptotic SD(MLE).
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
      n_param = sum(is.finite(z$calibration_rmse)),
      overall_calibration_rmse = sqrt(mean(z$calibration_rmse^2, na.rm = TRUE)),
      overall_efficiency_rmse_mle_scale = sqrt(mean(z$efficiency_rmse_mle_scale^2, na.rm = TRUE))
    )
  }))
  rownames(overall) <- NULL
  overall <- overall[order(overall$T, match(overall$method, names(method_cols))), ]

  pop_target <- data.frame(method = "MLE", population_target = 1, stringsAsFactors = FALSE)
  pop_target <- rbind(pop_target, data.frame(
    method = "GMM_Dtheta",
    population_target = sqrt(mean(pop$se_inflation_gmm_to_mle^2, na.rm = TRUE)),
    stringsAsFactors = FALSE
  ))
  if ("se_inflation_aug_to_mle" %in% names(pop) && any(is.finite(pop$se_inflation_aug_to_mle))) {
    pop_target <- rbind(pop_target, data.frame(
      method = "GMM_aug",
      population_target = sqrt(mean(pop$se_inflation_aug_to_mle^2, na.rm = TRUE)),
      stringsAsFactors = FALSE
    ))
  }

  list(parameter = d, overall = overall, population_target = pop_target, pop = pop)
}

plot_population_godambe_clear <- function(indir, popdir, outdir, prefix) {
  pop_path <- find_input("population_parameter_inflation.csv", indir, popdir, required = TRUE)
  eig_path <- find_input("population_relative_eigenvalues.csv", indir, popdir, required = TRUE)

  parinf <- read_csv_required(pop_path)
  eig <- read_csv_required(eig_path)

  need_pop <- c("param", "se_inflation_gmm_to_mle")
  need_eig <- c("eig_index", "covariance_eigenvalue")
  miss_pop <- setdiff(need_pop, names(parinf))
  miss_eig <- setdiff(need_eig, names(eig))
  if (length(miss_pop)) stopf("Missing columns in population_parameter_inflation.csv: ", paste(miss_pop, collapse = ", "))
  if (length(miss_eig)) stopf("Missing columns in population_relative_eigenvalues.csv: ", paste(miss_eig, collapse = ", "))

  short_method <- function(m) {
    out <- m
    out[m == "GMM_Dtheta"] <- "just-id Dtheta GMM"
    out[m == "GMM_aug"] <- "over-id bounded-inv GMM"
    out
  }

  parinf <- ensure_group(parinf)
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

  eig$eig_index <- as.numeric(eig$eig_index)
  eig$covariance_eigenvalue <- as.numeric(eig$covariance_eigenvalue)
  if (!"method" %in% names(eig)) eig$method <- "GMM_Dtheta"
  eig_methods <- intersect(c("GMM_Dtheta", "GMM_aug"), unique(eig$method))

  legend_methods <- unique(c(gmm_methods, eig_methods))
  legend_methods <- intersect(c("GMM_Dtheta", "GMM_aug"), legend_methods)

  out_png <- file.path(outdir, paste0(prefix, "_population_godambe_clear.png"))

  save_png(out_png, width = 1900, height = 1120, res = 170, draw_fun = function() {
    old <- par(no.readonly = TRUE)
    on.exit(par(old), add = TRUE)

    layout(
      matrix(c(1, 2,
               3, 3), nrow = 2, byrow = TRUE),
      widths = c(1.35, 1.0),
      heights = c(1.0, 0.13)
    )

    par(oma = c(0, 0, 3.0, 0))

    # --------------------------
    # Left panel: parameter-wise SE inflation
    # --------------------------
    par(mar = c(5.0, 7.8, 3.5, 1.2))

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
      xlab = "asymptotic SE ratio vs MLE",
      ylab = "",
      main = "Per-parameter SE inflation",
      cex.main = 1.25,
      cex.lab = 1.08,
      cex.axis = 0.98
    )

    axis(2, at = y, labels = params, las = 2, cex.axis = 0.98)
    grid(nx = NULL, ny = NA, col = "grey88")
    abline(v = 1, lty = 3, lwd = 1.8, col = "grey25")

    for (m in gmm_methods) {
      z <- par_long[par_long$method == m, , drop = FALSE]
      yy <- y[match(z$param, params)] + offsets[[m]]
      segments(
        1, yy, z$se_ratio, yy,
        col = method_cols[[m]],
        lwd = 5,
        lend = "butt"
      )
      points(
        z$se_ratio, yy,
        pch = method_pch[[m]],
        col = method_cols[[m]],
        cex = 1.1
      )
    }

    # --------------------------
    # Right panel: covariance eigenvalue inflation
    # --------------------------
    par(mar = c(5.0, 5.1, 3.5, 1.2))

    yr <- safe_range(eig$covariance_eigenvalue, include = 1, pad = 0.14)

    plot(
      NA,
      xlim = range(eig$eig_index),
      ylim = yr,
      xlab = "ordered eigenvalue index",
      ylab = "eigenvalue ratio vs MLE",
      main = "Covariance eigenvalue inflation",
      cex.main = 1.25,
      cex.lab = 1.08,
      cex.axis = 0.98
    )

    grid(col = "grey88")
    abline(h = 1, lty = 3, lwd = 1.8, col = "grey25")

    for (m in eig_methods) {
      z <- eig[eig$method == m, , drop = FALSE]
      z <- z[order(z$eig_index), ]
      lines(
        z$eig_index,
        z$covariance_eigenvalue,
        type = "b",
        pch = method_pch[[m]],
        lty = method_lty[[m]],
        lwd = 3,
        col = method_cols[[m]]
      )
    }

    # --------------------------
    # Shared legend below both panels
    # --------------------------
    par(mar = c(0, 0, 0, 0))
    plot.new()

    legend(
      "center",
      legend = c(short_method(legend_methods), "MLE benchmark"),
      col = c(method_cols[legend_methods], "grey25"),
      pch = c(method_pch[legend_methods], NA),
      lty = c(method_lty[legend_methods], 3),
      lwd = c(rep(3, length(legend_methods)), 1.8),
      bty = "n",
      cex = 0.90,
      ncol = length(legend_methods) + 1,
      x.intersp = 0.75,
      seg.len = 2.2
    )

    mtext(
      "Godambe efficiency loss relative to MLE",
      outer = TRUE,
      cex = 1.45,
      font = 2,
      line = 1.0
    )
  })

  out_png
}


plot_overall_efficiency <- function(overall, pop_target, outdir, prefix, tmin) {
  d <- overall[overall$T >= tmin, , drop = FALSE]
  if (!nrow(d)) stopf("No overall rows with T >= ", tmin)

  out_png <- file.path(outdir, paste0(prefix, "_overall_efficiency_Tge", tmin, ".png"))
  save_png(out_png, width = 1700, height = 1050, res = 170, draw_fun = function() {
    par(mar = c(5.0, 5.5, 4.9, 1.2))
    yy <- c(d$overall_efficiency_rmse_mle_scale, pop_target$population_target, 1)
    plot(NA, xlim = range(d$T), ylim = safe_range(yy, include = 1, pad = 0.12),
         log = "x", xaxt = "n",
         xlab = "T",
         ylab = "overall normalized RMSE",
         main = paste0("Overall normalized RMSE on common MLE scale, T >= ", tmin),
         cex.main = 1.35, cex.lab = 1.15, cex.axis = 1.05)
    mtext("sqrt(T) error normalized parameterwise by asymptotic SD(MLE); GMM target is the aggregate Godambe penalty.",
          side = 3, line = 0.55, cex = 0.92)
    log_x_axis(d$T)
    grid(col = "grey88")
    abline(h = 1, lty = 3, lwd = 1.7, col = "grey25")

    methods <- intersect(names(method_cols), unique(d$method))
    for (m in methods) {
      z <- d[d$method == m, , drop = FALSE]
      if (!nrow(z)) next
      z <- z[order(z$T), ]
      lines(z$T, z$overall_efficiency_rmse_mle_scale,
            type = "b",
            pch = method_pch[[m]],
            lty = method_lty[[m]],
            col = method_cols[[m]], lwd = 3.2, cex = 1.25)
      target <- pop_target$population_target[pop_target$method == m]
      if (length(target) && is.finite(target[1])) {
        abline(h = target[1], col = adjustcolor(method_cols[[m]], alpha.f = 0.45), lwd = 3)
      }
    }

    legend("topright",
           legend = c(paste0(pretty_method(methods), " empirical"),
                      "faint colored horizontal lines = population targets",
                      "black dotted horizontal line = parity"),
           col = c(method_cols[methods], "grey45", "grey25"),
           pch = c(method_pch[methods], NA, NA),
           lty = c(method_lty[methods], 1, 3),
           lwd = c(rep(3.2, length(methods)), 3, 1.7),
           bty = "n", cex = 0.86)
  })
  out_png
}

plot_overall_calibration <- function(overall, outdir, prefix, tmin) {
  d <- overall[overall$T >= tmin, , drop = FALSE]
  if (!nrow(d)) stopf("No overall rows with T >= ", tmin)

  out_png <- file.path(outdir, paste0(prefix, "_overall_calibration_Tge", tmin, ".png"))
  save_png(out_png, width = 1700, height = 1050, res = 170, draw_fun = function() {
    par(mar = c(5.0, 5.5, 4.9, 1.2))
    plot(NA, xlim = range(d$T), ylim = safe_range(d$overall_calibration_rmse, include = 1, pad = 0.12),
         log = "x", xaxt = "n",
         xlab = "T",
         ylab = "overall standardized RMSE",
         main = paste0("Overall covariance calibration, T >= ", tmin),
         cex.main = 1.35, cex.lab = 1.15, cex.axis = 1.05)
    mtext("sqrt(T) error normalized by each method's own asymptotic SD; target is 1 for both methods.",
          side = 3, line = 0.55, cex = 0.92)
    log_x_axis(d$T)
    grid(col = "grey88")
    abline(h = 1, lty = 3, lwd = 1.7, col = "grey25")

    methods <- intersect(names(method_cols), unique(d$method))
    for (m in methods) {
      z <- d[d$method == m, , drop = FALSE]
      if (!nrow(z)) next
      z <- z[order(z$T), ]
      lines(z$T, z$overall_calibration_rmse,
            type = "b",
            pch = method_pch[[m]],
            lty = method_lty[[m]],
            col = method_cols[[m]], lwd = 3.2, cex = 1.25)
    }

    legend("topright",
           legend = c(pretty_method(methods), "target = 1"),
           col = c(method_cols[methods], "grey25"),
           pch = c(method_pch[methods], NA),
           lty = c(method_lty[methods], 3),
           lwd = c(rep(3.2, length(methods)), 1.7),
           bty = "n", cex = 0.88)
  })
  out_png
}

plot_ci_coverage_width <- function(indir, popdir, outdir, prefix, tmin) {
  # Width-only plot using parameter-level CI summaries.
  path <- find_input("bivar_asymptotic_ci_by_parameter.csv", indir, popdir, required = FALSE)
  if (is.na(path)) {
    path <- find_input("bivar_asymptotic_ci_overall.csv", indir, popdir, required = FALSE)
  }
  if (is.na(path)) return(NA_character_)

  d <- read_csv_required(path)

  need <- c("T", "method")
  miss <- setdiff(need, names(d))
  if (length(miss)) {
    stopf("Missing columns in CI summary file: ", paste(miss, collapse = ", "))
  }

  d$T <- as.numeric(d$T)

  # Find an absolute CI-width column. If only half-width exists, double it.
  width_candidates <- c(
    "mean_width", "ci_width", "width", "avg_width",
    "mean_ci_width", "average_width", "full_width",
    "mean_full_width"
  )
  width_col <- intersect(width_candidates, names(d))

  if (length(width_col) >= 1) {
    d$plot_width <- as.numeric(d[[width_col[1]]])
  } else if ("half_width" %in% names(d)) {
    d$plot_width <- 2 * as.numeric(d$half_width)
  } else if ("mean_half_width" %in% names(d)) {
    d$plot_width <- 2 * as.numeric(d$mean_half_width)
  } else {
    stopf(
      "Could not find an absolute width column in ", basename(path),
      ". Available columns are: ", paste(names(d), collapse = ", ")
    )
  }

  d <- d[is.finite(d$T) & d$T >= tmin & is.finite(d$plot_width), , drop = FALSE]
  if (!nrow(d)) return(NA_character_)

  # Average over parameter coordinates if this is the by-parameter file.
  agg <- aggregate(plot_width ~ T + method, data = d, FUN = mean, na.rm = TRUE)
  names(agg)[names(agg) == "plot_width"] <- "mean_width"

  out_png <- file.path(outdir, paste0(prefix, "_ci_width_Tge", tmin, ".png"))

  save_png(out_png, width = 1400, height = 1150, res = 170, draw_fun = function() {
    old <- par(no.readonly = TRUE)
    on.exit(par(old), add = TRUE)

    par(mfrow = c(1, 1))
    par(oma = c(0, 0, 0, 0))
    par(mar = c(5.0, 5.4, 3.7, 1.4))

    methods <- intersect(names(method_cols), unique(agg$method))
    yy <- agg$mean_width
    ypad <- 0.06 * (max(yy, na.rm = TRUE) - min(yy, na.rm = TRUE) + 1e-8)

    plot(
      NA,
      xlim = range(agg$T),
      ylim = c(min(yy, na.rm = TRUE) - ypad, max(yy, na.rm = TRUE) + ypad),
      log = "x",
      xaxt = "n",
      xlab = "T",
      ylab = "mean asymptotic CI width",
      main = "Asymptotic 95% CI width",
      cex.main = 1.25,
      cex.lab = 1.05,
      cex.axis = 0.95
    )

    log_x_axis(agg$T)
    grid(col = "grey88")

    for (m in methods) {
      z <- agg[agg$method == m, , drop = FALSE]
      if (!nrow(z)) next
      z <- z[order(z$T), ]

      lines(
        z$T,
        z$mean_width,
        type = "b",
        pch = method_pch[[m]],
        lty = method_lty[[m]],
        col = method_cols[[m]],
        lwd = 3.0,
        cex = 1.15
      )
    }

    legend(
      "topright",
      legend = pretty_method(methods),
      col = method_cols[methods],
      pch = method_pch[methods],
      lty = method_lty[methods],
      lwd = rep(3, length(methods)),
      bty = "n",
      cex = 0.82
    )
  })

  out_png
}


plot_loglog_rmse <- function(indir, popdir, outdir, prefix) {
  path <- find_input("bivar_summary_by_parameter.csv", indir, popdir, required = TRUE)
  d <- read_csv_required(path)
  d <- ensure_group(d)

  need <- c("T", "method", "param", "group", "rmse")
  miss <- setdiff(need, names(d))
  if (length(miss)) stopf("Missing columns in bivar_summary_by_parameter.csv: ", paste(miss, collapse = ", "))

  d$T <- as.numeric(d$T)
  d$rmse <- as.numeric(d$rmse)
  d <- d[is.finite(d$T) & is.finite(d$rmse) & d$rmse > 0, , drop = FALSE]

  out_png <- file.path(outdir, paste0(prefix, "_loglog_rmse.png"))
  save_png(out_png, width = 1900, height = 1550, res = 170, draw_fun = function() {
    layout(matrix(c(1, 1, 2, 3, 4, 5), nrow = 3, byrow = TRUE), heights = c(0.18, 1, 1))

    par(mar = c(0, 0, 0, 0))
    plot.new()
    methods <- intersect(names(method_cols), unique(d$method))
    legend("center",
           legend = c(pretty_method(methods), "reference slope -1/2"),
           col = c(method_cols[methods], "grey25"),
           pch = c(method_pch[methods], NA),
           lty = c(method_lty[methods], 3),
           lwd = c(rep(3, length(methods)), 1.7),
           bty = "n", ncol = min(4, length(methods) + 1), cex = 0.95)

    title_map <- c(
      "baselines" = "Baselines",
      "self-excitation" = "Self-excitation",
      "cross-excitation" = "Cross-excitation",
      "decay" = "Decay"
    )

    for (g in group_order) {
      z <- d[d$group == g, , drop = FALSE]
      par(mar = c(4.7, 5.4, 3.2, 1.2))
      if (!nrow(z)) {
        plot.new()
        next
      }

      plot(NA,
           xlim = range(z$T),
           ylim = safe_range(z$rmse, pad = 0.10, lower_floor = .Machine$double.eps),
           log = "xy", xaxt = "n",
           xlab = "T", ylab = "RMSE",
           main = title_map[[g]],
           cex.main = 1.2, cex.lab = 1.05, cex.axis = 0.95)
      log_x_axis(z$T)
      grid(col = "grey88")

      for (m in methods) {
        for (p in intersect(param_order, unique(z$param))) {
          zz <- z[z$method == m & z$param == p, , drop = FALSE]
          if (!nrow(zz)) next
          zz <- zz[order(zz$T), ]
          lines(zz$T, zz$rmse,
                type = "b",
                col = method_cols[[m]],
                lty = method_lty[[m]],
                pch = param_pch[[p]],
                lwd = 2.4, cex = 1.0)
        }
      }

      TT <- sort(unique(z$T))
      y0 <- median(z$rmse[z$T == min(TT)], na.rm = TRUE)
      lines(TT, y0 * (TT / min(TT))^(-0.5), lty = 3, lwd = 1.7, col = "grey25")

      ps <- intersect(param_order, unique(z$param))
      legend("topright",
             legend = ps,
             col = "grey25",
             pch = param_pch[ps],
             lty = param_lty[ps],
             bty = "n", cex = 0.76)
    }

    mtext("Log-log RMSE curves", outer = TRUE, line = -1.2, cex = 1.45, font = 2)
  })
  out_png
}

main <- function() {
  opts <- parse_args()

  indir <- normalizePath(opts$indir, mustWork = FALSE)
  popdir <- opts$popdir
  if (is.na(popdir) || !nzchar(popdir)) popdir <- NA_character_ else popdir <- normalizePath(popdir, mustWork = FALSE)

  outdir <- opts$outdir
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  pack <- make_overall_tables(indir, popdir)

  outputs <- c(
    population_godambe = plot_population_godambe_clear(indir, popdir, outdir, opts$prefix),
    overall_efficiency = plot_overall_efficiency(pack$overall, pack$population_target, outdir, opts$prefix, opts$tmin),
    overall_calibration = plot_overall_calibration(pack$overall, outdir, opts$prefix, opts$tmin),
    ci_coverage_width = plot_ci_coverage_width(indir, popdir, outdir, opts$prefix, opts$tmin),
    loglog_rmse = plot_loglog_rmse(indir, popdir, outdir, opts$prefix)
  )

  manifest <- data.frame(
    figure = names(outputs),
    file = unname(outputs),
    stringsAsFactors = FALSE
  )
  write.csv(manifest, file.path(outdir, paste0(opts$prefix, "_plot_manifest.csv")), row.names = FALSE)
  # Backward-compatible manifest name used by earlier scripts.
  write.csv(manifest, file.path(outdir, paste0(opts$prefix, "_four_plot_manifest.csv")), row.names = FALSE)

  cat("Wrote paper plots to:\n")
  cat("  ", normalizePath(outdir, mustWork = FALSE), "\n\n", sep = "")
  print(manifest, row.names = FALSE)
}

main()

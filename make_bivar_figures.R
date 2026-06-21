#!/usr/bin/env Rscript
################################################################################
# Base-R figure maker for the focused bivariate Hawkes Godambe study.
# Uses a colorblind-safe palette and explicit legends.
################################################################################

parse_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list(); i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (!startsWith(a, "--")) { i <- i + 1L; next }
    a <- substring(a, 3L)
    if (grepl("=", a, fixed = TRUE)) { out[[sub("=.*$", "", a)]] <- sub("^[^=]*=", "", a); i <- i + 1L }
    else { key <- a; if (i < length(args) && !startsWith(args[[i + 1L]], "--")) { out[[key]] <- args[[i + 1L]]; i <- i + 2L } else { out[[key]] <- "TRUE"; i <- i + 1L } }
  }
  out
}
arg <- function(opts, name, default = NULL) if (!is.null(opts[[name]])) opts[[name]] else default
read_or_null <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0) return(NULL)
  tryCatch({
    x <- read.csv(path, stringsAsFactors = FALSE)
    if (!nrow(x) && !ncol(x)) NULL else x
  }, error = function(e) NULL)
}

param_order <- c("mu1", "mu2", "alpha11", "alpha12", "alpha21", "alpha22", "beta")
group_order <- c("baselines", "self-excitation", "cross-excitation", "decay")
param_group <- function(p) {
  out <- rep("other", length(p)); names(out) <- p
  out[p %in% c("mu1", "mu2")] <- "baselines"
  out[p %in% c("alpha11", "alpha22")] <- "self-excitation"
  out[p %in% c("alpha12", "alpha21")] <- "cross-excitation"
  out[p %in% "beta"] <- "decay"
  out
}
pretty_method <- function(x) {
  out <- x
  out[x == "GMM_Dtheta"] <- "just-identified Dtheta GMM"
  out[x == "GMM_aug"] <- "overidentified bounded-inverse GMM"
  out
}
method_colors <- c(MLE = "#0072B2", GMM_Dtheta = "#D55E00", GMM_aug = "#009E73")
method_lty <- c(MLE = 1, GMM_Dtheta = 2, GMM_aug = 4)
method_pch <- c(MLE = 16, GMM_Dtheta = 17, GMM_aug = 15)
group_colors <- c("baselines" = "#009E73", "self-excitation" = "#CC79A7", "cross-excitation" = "#E69F00", "decay" = "#56B4E9")
prob_pch <- c("0.5" = 16, "0.9" = 17, "0.95" = 15, "0.99" = 18)
param_lty <- setNames(c(1, 2, 1, 2, 1, 2, 1), param_order)
param_pch <- setNames(c(16, 17, 16, 17, 16, 17, 15), param_order)

safe_ylim <- function(y, include = NULL, pad = 0.08, positive = FALSE) {
  y <- c(y, include); y <- y[is.finite(y)]
  if (!length(y)) return(if (positive) c(0.1, 1) else c(0, 1))
  r <- range(y)
  if (positive) r[1] <- max(min(r[1], r[2] * 0.8), .Machine$double.eps)
  if (diff(r) == 0) r <- r + c(-0.5, 0.5) * max(1, abs(r[1]))
  r + c(-pad, pad) * diff(r)
}
add_xaxis <- function(T) { ticks <- sort(unique(as.numeric(T))); axis(1, at = ticks, labels = ticks) }

plot_population_godambe <- function(indir, outdir) {
  eig <- read_or_null(file.path(indir, "population_relative_eigenvalues.csv"))
  summ <- read_or_null(file.path(indir, "population_godambe_summary.csv"))
  if (is.null(eig) || !nrow(eig)) return(FALSE)
  pdf(file.path(outdir, "population_godambe_eigenvalues.pdf"), width = 8.2, height = 4.8)
  old <- par(no.readonly = TRUE); on.exit({ par(old); dev.off() }, add = TRUE)
  layout(matrix(c(1, 2), nrow = 1), widths = c(1.35, 1))
  par(mar = c(4.4, 5.2, 2.4, 1.0))
  if (!"method" %in% names(eig)) eig$method <- "GMM_Dtheta"
  plot(NA, xlim = range(eig$eig_index), ylim = safe_ylim(eig$covariance_eigenvalue, include = 1),
       xlab = "ordered eigenvalue index",
       ylab = expression("eigenvalues of "*I^{1/2}*V*I^{1/2}),
       main = "Population Godambe inflation")
  abline(h = 1, lty = 3, lwd = 1.5)
  grid(col = "grey85")
  for (m in intersect(names(method_colors), unique(eig$method))) {
    zz <- eig[eig$method == m, ]; zz <- zz[order(zz$eig_index), ]
    lines(zz$eig_index, zz$covariance_eigenvalue, type = "b", pch = method_pch[[m]],
          lty = method_lty[[m]], lwd = 2, col = method_colors[[m]])
  }
  legend("topleft", legend = c(pretty_method(intersect(names(method_colors), unique(eig$method))), "MLE benchmark"),
         col = c(method_colors[intersect(names(method_colors), unique(eig$method))], "black"),
         pch = c(method_pch[intersect(names(method_colors), unique(eig$method))], NA),
         lty = c(method_lty[intersect(names(method_colors), unique(eig$method))], 3),
         lwd = c(rep(2, length(intersect(names(method_colors), unique(eig$method)))), 1.5), bty = "n", cex = 0.82)

  par(mar = c(6.8, 5.0, 2.4, 1.0))
  if (!is.null(summ) && nrow(summ)) {
    stats <- c("eig_min_ls_to_mle", "eig_median_ls_to_mle", "eig_max_ls_to_mle", "trace_ratio_ls_to_mle",
               "eig_min_aug_to_mle", "eig_median_aug_to_mle", "eig_max_aug_to_mle", "trace_ratio_aug_to_mle")
    labs <- c("D min eig", "D median", "D max", "D trace", "Aug min", "Aug median", "Aug max", "Aug trace")
    dd <- summ[summ$statistic %in% stats, ]
    dd <- dd[match(stats, dd$statistic), ]
    xx <- seq_len(nrow(dd))
    plot(NA, xlim = c(0.5, length(xx) + 0.5), ylim = safe_ylim(c(dd$estimate, dd$lo, dd$hi), include = 1),
         xaxt = "n", xlab = "", ylab = "relative to MLE")
    abline(h = 1, lty = 3, lwd = 1.5)
    grid(col = "grey88")
    axis(1, at = xx, labels = labs, las = 2)
    if (all(c("lo", "hi") %in% names(dd)) && any(is.finite(dd$lo))) {
      arrows(xx, dd$lo, xx, dd$hi, angle = 90, code = 3, length = 0.04, col = "grey35")
    }
    points(xx, dd$estimate, pch = 16, cex = 1.1, col = "#D55E00")
    title("Summary ratios")
  } else plot.new()
  TRUE
}

plot_parameter_inflation <- function(indir, outdir) {
  dat_long <- read_or_null(file.path(indir, "population_parameter_inflation_long.csv"))
  if (!is.null(dat_long) && nrow(dat_long)) {
    dat <- dat_long
  } else {
    wide <- read_or_null(file.path(indir, "population_parameter_inflation.csv")); if (is.null(wide) || !nrow(wide)) return(FALSE)
    dat <- data.frame(method = "GMM_Dtheta", param = wide$param, group = wide$group,
                      se_inflation_to_mle = wide$se_inflation_gmm_to_mle, stringsAsFactors = FALSE)
    if ("se_inflation_aug_to_mle" %in% names(wide)) {
      dat <- rbind(dat, data.frame(method = "GMM_aug", param = wide$param, group = wide$group,
                                   se_inflation_to_mle = wide$se_inflation_aug_to_mle, stringsAsFactors = FALSE))
    }
  }
  dat <- dat[is.finite(dat$se_inflation_to_mle), ]
  if (!nrow(dat)) return(FALSE)
  dat$param <- factor(dat$param, levels = param_order)
  dat <- dat[order(dat$param, dat$method), ]
  methods <- intersect(names(method_colors), unique(dat$method))
  mat <- sapply(methods, function(m) {
    out <- rep(NA_real_, length(param_order)); names(out) <- param_order
    z <- dat[dat$method == m, ]
    out[as.character(z$param)] <- z$se_inflation_to_mle
    out
  })
  pdf(file.path(outdir, "asymptotic_se_inflation_by_parameter.pdf"), width = 8.4, height = 4.8)
  old <- par(no.readonly = TRUE); on.exit({ par(old); dev.off() }, add = TRUE)
  par(mar = c(5.2, 5.0, 2.2, 1.0))
  bp <- barplot(t(mat), beside = TRUE, names.arg = param_order,
                col = method_colors[methods], border = "grey30",
                ylim = safe_ylim(mat, include = 1),
                ylab = "asymptotic SE ratio: GMM / MLE", main = "Per-parameter Godambe penalty")
  abline(h = 1, lty = 3, lwd = 1.5)
  grid(nx = NA, ny = NULL, col = "grey88")
  legend("topright", legend = pretty_method(methods), fill = method_colors[methods], bty = "n", cex = 0.82)
  TRUE
}

plot_standardized_rmse <- function(indir, outdir) {
  dat <- read_or_null(file.path(indir, "bivar_standardized_summary.csv")); if (is.null(dat) || !nrow(dat)) return(FALSE)
  dat <- dat[is.finite(dat$scaled_rmse_over_target), ]
  if (!nrow(dat)) return(FALSE)
  pdf(file.path(outdir, "standardized_rmse_by_parameter.pdf"), width = 8.5, height = 8.0)
  old <- par(no.readonly = TRUE); on.exit({ par(old); dev.off() }, add = TRUE)
  layout(matrix(c(1, 1, 2, 3, 4, 5), ncol = 2, byrow = TRUE), heights = c(0.35, 1, 1))
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", legend = c(pretty_method(names(method_colors)), "target = 1"),
         col = c(method_colors, "black"), pch = c(method_pch, NA), lty = c(method_lty, 3), lwd = c(rep(2, length(method_colors)), 1.5), bty = "n", ncol = 3)
  for (g in group_order) {
    ss <- dat[dat$group == g, ]
    par(mar = c(4.3, 5.1, 2.2, 1.0))
    plot(NA, xlim = range(ss$T), ylim = safe_ylim(ss$scaled_rmse_over_target, include = 1), log = "x", xaxt = "n",
         xlab = "T", ylab = expression(sqrt(T)*" RMSE / asymptotic SD"), main = g)
    add_xaxis(ss$T); abline(h = 1, lty = 3, lwd = 1.5); grid(col = "grey88")
    for (m in intersect(names(method_colors), unique(ss$method))) {
      for (p in intersect(param_order, unique(ss$param))) {
        mm <- ss[ss$method == m & ss$param == p, ]; if (!nrow(mm)) next
        mm <- mm[order(mm$T), ]
        lines(mm$T, mm$scaled_rmse_over_target, type = "b", pch = param_pch[[p]], lty = param_lty[[p]],
              col = method_colors[[m]], lwd = 2)
      }
    }
    legend("topright", legend = unique(as.character(ss$param)), pch = param_pch[unique(as.character(ss$param))],
           lty = param_lty[unique(as.character(ss$param))], col = "grey25", bty = "n", cex = 0.82)
  }
  TRUE
}

plot_loglog_rmse <- function(indir, outdir) {
  dat <- read_or_null(file.path(indir, "bivar_summary_by_parameter.csv")); if (is.null(dat) || !nrow(dat)) return(FALSE)
  dat <- dat[dat$rmse > 0 & is.finite(dat$rmse), ]
  if (!nrow(dat)) return(FALSE)
  pdf(file.path(outdir, "loglog_rmse_curves.pdf"), width = 8.5, height = 8.0)
  old <- par(no.readonly = TRUE); on.exit({ par(old); dev.off() }, add = TRUE)
  layout(matrix(c(1, 1, 2, 3, 4, 5), ncol = 2, byrow = TRUE), heights = c(0.35, 1, 1))
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", legend = c(pretty_method(names(method_colors)), "reference slope -1/2"),
         col = c(method_colors, "black"), pch = c(method_pch, NA), lty = c(method_lty, 3), lwd = c(rep(2, length(method_colors)), 1.5), bty = "n", ncol = 3)
  for (g in group_order) {
    ss <- dat[dat$group == g, ]
    par(mar = c(4.3, 5.1, 2.2, 1.0))
    plot(NA, xlim = range(ss$T), ylim = safe_ylim(ss$rmse, positive = TRUE), log = "xy", xaxt = "n",
         xlab = "T", ylab = "RMSE", main = g)
    add_xaxis(ss$T); grid(col = "grey88")
    for (m in intersect(names(method_colors), unique(ss$method))) {
      for (p in intersect(param_order, unique(ss$param))) {
        mm <- ss[ss$method == m & ss$param == p, ]; if (!nrow(mm)) next
        mm <- mm[order(mm$T), ]
        lines(mm$T, mm$rmse, type = "b", pch = param_pch[[p]], lty = param_lty[[p]], col = method_colors[[m]], lwd = 2)
      }
    }
    # Reference slope anchored at the median first-T RMSE in the panel.
    TT <- sort(unique(ss$T)); y0 <- median(ss$rmse[ss$T == min(TT)], na.rm = TRUE)
    lines(TT, y0 * (TT / min(TT))^(-0.5), lty = 3, lwd = 1.5)
    legend("topright", legend = unique(as.character(ss$param)), pch = param_pch[unique(as.character(ss$param))],
           lty = param_lty[unique(as.character(ss$param))], col = "grey25", bty = "n", cex = 0.82)
  }
  TRUE
}

plot_rate_slopes <- function(indir, outdir) {
  dat <- read_or_null(file.path(indir, "bivar_rate_slopes.csv")); if (is.null(dat) || !nrow(dat)) return(FALSE)
  dat$param <- factor(dat$param, levels = param_order)
  dat <- dat[order(dat$group, dat$param, dat$method), ]
  pdf(file.path(outdir, "loglog_rmse_slopes.pdf"), width = 8.8, height = 5.5)
  old <- par(no.readonly = TRUE); on.exit({ par(old); dev.off() }, add = TRUE)
  par(mar = c(7.2, 5.0, 2.2, 1.0))
  keys <- paste(pretty_method(dat$method), dat$param, sep = ": ")
  xx <- seq_len(nrow(dat))
  cols <- method_colors[dat$method]
  plot(xx, dat$slope_log_rmse_on_log_T, xaxt = "n", xlab = "", ylab = "slope of log RMSE on log T",
       ylim = safe_ylim(dat$slope_log_rmse_on_log_T, include = -0.5), pch = 16, col = cols, cex = 1.1,
       main = "Rate check: target slope -1/2")
  abline(h = -0.5, lty = 3, lwd = 1.5)
  grid(nx = NA, ny = NULL, col = "grey88")
  axis(1, at = xx, labels = keys, las = 2, cex.axis = 0.72)
  legend("topright", legend = pretty_method(names(method_colors)), col = method_colors, pch = 16, bty = "n")
  TRUE
}

plot_lan_quantiles <- function(indir, outdir) {
  dat <- read_or_null(file.path(indir, "bivar_lan_quadratic_quantiles.csv")); if (is.null(dat) || !nrow(dat)) return(FALSE)
  dat <- dat[dat$prob %in% c(0.90, 0.95, 0.99), ]; if (!nrow(dat)) return(FALSE)
  pdf(file.path(outdir, "lan_quadratic_quantile_ratios.pdf"), width = 8.0, height = 5.2)
  old <- par(no.readonly = TRUE); on.exit({ par(old); dev.off() }, add = TRUE)
  par(mar = c(4.4, 5.1, 2.2, 1.0))
  plot(NA, xlim = range(dat$T), ylim = safe_ylim(dat$ratio_to_chisq, include = 1), log = "x", xaxt = "n",
       xlab = "T", ylab = "empirical / chi-square quantile", main = "LAN/Godambe quadratic-form calibration")
  add_xaxis(dat$T); abline(h = 1, lty = 3, lwd = 1.5); grid(col = "grey88")
  leg <- character(); leg_col <- character(); leg_pch <- integer(); leg_lty <- integer()
  for (m in intersect(names(method_colors), unique(dat$method))) {
    for (pr in sort(unique(dat$prob))) {
      mm <- dat[dat$method == m & abs(dat$prob - pr) < 1e-12, ]; if (!nrow(mm)) next
      mm <- mm[order(mm$T), ]
      pch <- prob_pch[[as.character(pr)]]; if (is.null(pch)) pch <- 16
      lines(mm$T, mm$ratio_to_chisq, type = "b", pch = pch, lty = method_lty[[m]],
            col = method_colors[[m]], lwd = 2)
      leg <- c(leg, paste(pretty_method(m), "q=", pr)); leg_col <- c(leg_col, method_colors[[m]]); leg_pch <- c(leg_pch, pch); leg_lty <- c(leg_lty, method_lty[[m]])
    }
  }
  legend("topright", legend = leg, col = leg_col, pch = leg_pch, lty = leg_lty, lwd = 2, bty = "n", cex = 0.78)
  TRUE
}

plot_asymptotic_ci_coverage_width <- function(indir, outdir) {
  overall <- read_or_null(file.path(indir, "bivar_asymptotic_ci_overall.csv"))
  if (is.null(overall) || !nrow(overall)) return(FALSE)

  overall$T <- as.numeric(overall$T)
  overall$mean_width <- as.numeric(overall$mean_width)
  overall <- overall[is.finite(overall$T) & is.finite(overall$mean_width), , drop = FALSE]
  if (!nrow(overall)) return(FALSE)

  pdf(file.path(outdir, "asymptotic_ci_width.pdf"), width = 7.2, height = 5.8)
  old <- par(no.readonly = TRUE)
  on.exit({ par(old); dev.off() }, add = TRUE)

  par(mfrow = c(1, 1))
  par(mar = c(4.4, 5.1, 2.6, 1.0))

  methods <- intersect(names(method_colors), unique(overall$method))
  yy <- overall$mean_width
  ypad <- 0.06 * (max(yy, na.rm = TRUE) - min(yy, na.rm = TRUE) + 1e-8)

  plot(
    NA,
    xlim = range(overall$T),
    ylim = c(min(yy, na.rm = TRUE) - ypad, max(yy, na.rm = TRUE) + ypad),
    log = "x",
    xaxt = "n",
    xlab = "T",
    ylab = "mean asymptotic CI width",
    main = "Asymptotic 95% CI width"
  )

  add_xaxis(overall$T)
  grid(col = "grey88")

  for (m in methods) {
    z <- overall[overall$method == m, , drop = FALSE]
    if (!nrow(z)) next
    z <- z[order(z$T), ]
    lines(
      z$T,
      z$mean_width,
      type = "b",
      pch = method_pch[[m]],
      lty = method_lty[[m]],
      col = method_colors[[m]],
      lwd = 2
    )
  }

  legend(
    "topright",
    legend = pretty_method(methods),
    col = method_colors[methods],
    pch = method_pch[methods],
    lty = method_lty[methods],
    lwd = rep(2, length(methods)),
    bty = "n",
    cex = 0.82
  )

  TRUE
}


plot_emp_cov <- function(indir, outdir) {
  dat <- read_or_null(file.path(indir, "bivar_empirical_vs_population_covariance.csv")); if (is.null(dat) || !nrow(dat)) return(FALSE)
  dat <- dat[dat$statistic %in% c("trace_ratio", "eig_median"), ]; if (!nrow(dat)) return(FALSE)
  stat_cols <- c(trace_ratio = "#009E73", eig_median = "#CC79A7")
  pdf(file.path(outdir, "empirical_vs_population_covariance.pdf"), width = 8.0, height = 5.2)
  old <- par(no.readonly = TRUE); on.exit({ par(old); dev.off() }, add = TRUE)
  par(mar = c(4.4, 5.1, 2.2, 1.0))
  plot(NA, xlim = range(dat$T), ylim = safe_ylim(c(dat$estimate, dat$lo, dat$hi), include = 1), log = "x", xaxt = "n",
       xlab = "T", ylab = "empirical covariance / population target", main = "Empirical covariance calibration")
  add_xaxis(dat$T); abline(h = 1, lty = 3, lwd = 1.5); grid(col = "grey88")
  for (m in intersect(names(method_colors), unique(dat$method))) {
    for (st in names(stat_cols)) {
      mm <- dat[dat$method == m & dat$statistic == st, ]; if (!nrow(mm)) next
      mm <- mm[order(mm$T), ]
      lines(mm$T, mm$estimate, type = "b", pch = ifelse(st == "trace_ratio", 16, 17),
            lty = method_lty[[m]], col = stat_cols[[st]], lwd = 2)
      if (all(c("lo", "hi") %in% names(mm)) && any(is.finite(mm$lo))) {
        arrows(mm$T, mm$lo, mm$T, mm$hi, angle = 90, code = 3, length = 0.03, col = adjustcolor(stat_cols[[st]], alpha.f = 0.5))
      }
    }
  }
  legend("topright", legend = c("trace ratio", "median eigenvalue", pretty_method(names(method_colors))),
         col = c(stat_cols, method_colors), pch = c(16, 17, rep(NA, length(method_colors))),
         lty = c(NA, NA, method_lty), lwd = c(NA, NA, rep(2, length(method_colors))), bty = "n", cex = 0.82)
  TRUE
}

main <- function() {
  opts <- parse_cli(); indir <- arg(opts, "indir", "bivar_results"); outdir <- arg(opts, "outdir", "bivar_figures")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  made <- c(
    population_godambe = plot_population_godambe(indir, outdir),
    parameter_inflation = plot_parameter_inflation(indir, outdir),
    standardized_rmse = plot_standardized_rmse(indir, outdir),
    loglog_rmse = plot_loglog_rmse(indir, outdir),
    rate_slopes = plot_rate_slopes(indir, outdir),
    lan_quantiles = plot_lan_quantiles(indir, outdir),
    ci_coverage_width = plot_asymptotic_ci_coverage_width(indir, outdir),
    empirical_covariance = plot_emp_cov(indir, outdir)
  )
  write.csv(data.frame(figure = names(made), created = as.logical(made)), file.path(outdir, "figure_manifest.csv"), row.names = FALSE)
}
main()

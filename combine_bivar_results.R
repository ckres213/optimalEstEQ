#!/usr/bin/env Rscript
################################################################################
# Combine focused bivariate Hawkes outputs and create analysis CSVs.
################################################################################

parse_cli0 <- function(args = commandArgs(trailingOnly = TRUE)) {
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
script_dir0 <- function() { ff <- commandArgs(FALSE); ff <- ff[grepl("^--file=", ff)]; if (length(ff)) dirname(normalizePath(sub("^--file=", "", ff[1L]))) else getwd() }
arg0 <- function(opts, name, default = NULL) if (!is.null(opts[[name]])) opts[[name]] else default
parse_int0 <- function(x, default) if (is.null(x) || !nzchar(as.character(x))) default else as.integer(x)
parse_num0 <- function(x, default) if (is.null(x) || !nzchar(as.character(x))) default else as.numeric(x)
parse_bool0 <- function(x, default = FALSE) { if (is.null(x) || !nzchar(as.character(x))) return(default); tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y", "on") }

load_core <- function(base_dir, tmp_build = TRUE) {
  core <- file.path(base_dir, "bivar_core.R"); cpp <- file.path(base_dir, "hawkes_bivar_fast.cpp")
  if (!file.exists(core) || !file.exists(cpp)) stop("Missing bivar_core.R or hawkes_bivar_fast.cpp in ", base_dir)
  old <- getwd()
  if (tmp_build) {
    tmp <- Sys.getenv("TMPDIR", unset = tempdir())
    bdir <- file.path(tmp, paste0("bivar_combine_build_", Sys.getpid(), "_", sample.int(1e9, 1)))
    dir.create(bdir, recursive = TRUE, showWarnings = FALSE)
    file.copy(core, file.path(bdir, "bivar_core.R"), overwrite = TRUE)
    file.copy(cpp, file.path(bdir, "hawkes_bivar_fast.cpp"), overwrite = TRUE)
    setwd(bdir); on.exit(setwd(old), add = TRUE); source("bivar_core.R")
  } else { setwd(base_dir); on.exit(setwd(old), add = TRUE); source("bivar_core.R") }
  invisible(TRUE)
}

read_csvs <- function(files) {
  files <- files[file.exists(files) & file.info(files)$size > 0]
  if (!length(files)) return(data.frame())
  parts <- lapply(files, function(f) tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL))
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts)) data.frame() else do.call(rbind, parts)
}

main <- function() {
  opts <- parse_cli0()
  base_dir <- normalizePath(arg0(opts, "base-dir", script_dir0()))
  indir <- arg0(opts, "indir", "bivar_results")
  popdir <- arg0(opts, "popdir", "bivar_population")
  outdir <- arg0(opts, "outdir", indir)
  B_boot <- parse_int0(arg0(opts, "B-boot", "500"), 500L)
  ci_level <- parse_num0(arg0(opts, "ci-level", "0.95"), 0.95)
  seed <- parse_int0(arg0(opts, "seed", "20260525"), 20260525L)
  tmp_build <- parse_bool0(arg0(opts, "tmp-build", "TRUE"), TRUE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  load_core(base_dir, tmp_build)

  spec <- bivar_spec()
  design <- data.frame(
    scenario = spec$name, description = spec$description, A = spec$A,
    mu1 = spec$theta[["mu1"]], mu2 = spec$theta[["mu2"]],
    alpha11 = spec$theta[["alpha11"]], alpha12 = spec$theta[["alpha12"]],
    alpha21 = spec$theta[["alpha21"]], alpha22 = spec$theta[["alpha22"]],
    beta = spec$theta[["beta"]], rho = spec$rho,
    lambda_bar1 = spec$lambda_bar[1], lambda_bar2 = spec$lambda_bar[2],
    burnin = burnin_length(spec$A, spec$rho), stringsAsFactors = FALSE)
  write.csv(design, file.path(outdir, "bivar_design.csv"), row.names = FALSE)

  pop_summary_files <- list.files(popdir, pattern = "^population_summary_array.*\\.csv$", full.names = TRUE)
  pop_matrix_files <- list.files(popdir, pattern = "^population_matrices_array.*\\.csv$", full.names = TRUE)
  pop <- NULL; pop_mat <- data.frame()
  if (length(pop_matrix_files)) {
    pop_mat <- read_csvs(pop_matrix_files)
    if (nrow(pop_mat)) {
      pop <- build_population_list(pop_mat)
      saveRDS(pop, file.path(outdir, "population_matrices.rds"))
      write.csv(population_summary_table(pop, B_boot = B_boot, seed = seed, mat_long = pop_mat),
                file.path(outdir, "population_godambe_summary.csv"), row.names = FALSE)
      write.csv(population_parameter_table(pop), file.path(outdir, "population_parameter_inflation.csv"), row.names = FALSE)
      write.csv(population_parameter_long_table(pop), file.path(outdir, "population_parameter_inflation_long.csv"), row.names = FALSE)
      write.csv(population_eigen_table(pop), file.path(outdir, "population_relative_eigenvalues.csv"), row.names = FALSE)
      if (length(pop_summary_files)) {
        ps <- read_csvs(pop_summary_files)
        if (nrow(ps)) write.csv(ps, file.path(outdir, "population_raw_summaries.csv"), row.names = FALSE)
      }
    }
  } else if (file.exists(file.path(outdir, "population_matrices.rds"))) {
    pop <- readRDS(file.path(outdir, "population_matrices.rds"))
    write.csv(population_parameter_table(pop), file.path(outdir, "population_parameter_inflation.csv"), row.names = FALSE)
    write.csv(population_parameter_long_table(pop), file.path(outdir, "population_parameter_inflation_long.csv"), row.names = FALSE)
    write.csv(population_eigen_table(pop), file.path(outdir, "population_relative_eigenvalues.csv"), row.names = FALSE)
  }

  raw_files <- list.files(indir, pattern = "^raw_array.*\\.csv$", full.names = TRUE)
  res <- read_csvs(raw_files)
  if (!nrow(res)) stop("No estimator raw_array*.csv files found in ", indir)
  key <- paste(res$T, res$rep, res$method, res$param, sep = "|")
  res <- res[!duplicated(key, fromLast = TRUE), , drop = FALSE]
  rownames(res) <- NULL
  write.csv(res, file.path(outdir, "bivar_raw_combined.csv"), row.names = FALSE)
  write.csv(summary_by_parameter(res), file.path(outdir, "bivar_summary_by_parameter.csv"), row.names = FALSE)
  write.csv(rate_slope_table(res), file.path(outdir, "bivar_rate_slopes.csv"), row.names = FALSE)

  if (!is.null(pop)) {
    write.csv(standardized_summary(res, pop), file.path(outdir, "bivar_standardized_summary.csv"), row.names = FALSE)
    write.csv(asymptotic_ci_long(res, pop, level = ci_level),
              file.path(outdir, "bivar_asymptotic_ci_long.csv"), row.names = FALSE)
    write.csv(asymptotic_ci_by_parameter(res, pop, level = ci_level),
              file.path(outdir, "bivar_asymptotic_ci_by_parameter.csv"), row.names = FALSE)
    write.csv(asymptotic_ci_overall(res, pop, level = ci_level),
              file.path(outdir, "bivar_asymptotic_ci_overall.csv"), row.names = FALSE)
    write.csv(lan_quadratic_quantiles(res, pop), file.path(outdir, "bivar_lan_quadratic_quantiles.csv"), row.names = FALSE)
    write.csv(empirical_vs_population_covariance(res, pop, B_boot = B_boot, seed = seed),
              file.path(outdir, "bivar_empirical_vs_population_covariance.csv"), row.names = FALSE)
  } else {
    message("Population matrices not found; skipped standardized/Godambe diagnostics")
  }

  fail_files <- c(list.files(indir, pattern = "^failures_array.*\\.csv$", full.names = TRUE),
                  list.files(popdir, pattern = "^population_failures_array.*\\.csv$", full.names = TRUE))
  fails <- read_csvs(fail_files)
  if (nrow(fails)) write.csv(fails, file.path(outdir, "bivar_failures_combined.csv"), row.names = FALSE)
  invisible(NULL)
}

main()

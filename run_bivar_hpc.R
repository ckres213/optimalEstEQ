#!/usr/bin/env Rscript
################################################################################
# Array runner for the focused bivariate Hawkes Godambe study.
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
arg0 <- function(opts, name, default = NULL, env = character()) { if (!is.null(opts[[name]])) return(opts[[name]]); for (e in env) { v <- Sys.getenv(e, unset = ""); if (nzchar(v)) return(v) }; default }
parse_int0 <- function(x, default) if (is.null(x) || !nzchar(as.character(x))) default else as.integer(x)
parse_num0 <- function(x, default) if (is.null(x) || !nzchar(as.character(x))) default else as.numeric(x)
parse_num_vec0 <- function(x, default) { if (is.null(x) || !nzchar(as.character(x))) return(default); z <- suppressWarnings(as.numeric(strsplit(as.character(x), ",")[[1L]])); z <- z[is.finite(z)]; if (!length(z)) default else z }
parse_bool0 <- function(x, default = FALSE) { if (is.null(x) || !nzchar(as.character(x))) return(default); tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y", "on") }

load_core <- function(base_dir, tmp_build = TRUE) {
  core <- file.path(base_dir, "bivar_core.R"); cpp <- file.path(base_dir, "hawkes_bivar_fast.cpp")
  if (!file.exists(core) || !file.exists(cpp)) stop("Missing bivar_core.R or hawkes_bivar_fast.cpp in ", base_dir)
  old <- getwd()
  if (tmp_build) {
    tmp <- Sys.getenv("TMPDIR", unset = tempdir())
    bdir <- file.path(tmp, paste0("bivar_hawkes_build_", Sys.getpid(), "_", sample.int(1e9, 1)))
    dir.create(bdir, recursive = TRUE, showWarnings = FALSE)
    file.copy(core, file.path(bdir, "bivar_core.R"), overwrite = TRUE)
    file.copy(cpp, file.path(bdir, "hawkes_bivar_fast.cpp"), overwrite = TRUE)
    setwd(bdir); on.exit(setwd(old), add = TRUE); source("bivar_core.R")
  } else { setwd(base_dir); on.exit(setwd(old), add = TRUE); source("bivar_core.R") }
  invisible(TRUE)
}

main <- function() {
  opts <- parse_cli0()
  base_dir <- normalizePath(arg0(opts, "base-dir", script_dir0()))
  outdir <- arg0(opts, "outdir", "bivar_results")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  array_id <- parse_int0(arg0(opts, "array-id", NULL, c("SLURM_ARRAY_TASK_ID", "PBS_ARRAY_INDEX", "PBS_ARRAYID")), 1L)
  array_count <- parse_int0(arg0(opts, "array-count", NULL, c("SLURM_ARRAY_TASK_COUNT")), 1L)
  if (is.na(array_count) || array_count < 1L) array_count <- parse_int0(arg0(opts, "n-array", "1"), 1L)
  workers <- parse_int0(arg0(opts, "workers", NULL, c("SLURM_CPUS_PER_TASK", "PBS_NCPUS")), 1L)
  T_grid <- parse_num_vec0(arg0(opts, "T-grid", "500,1000,2000,4000,8000"), c(500, 1000, 2000, 4000, 8000))
  reps_per_T <- parse_int0(arg0(opts, "reps-per-T", "400"), 400L)
  seed <- parse_int0(arg0(opts, "seed", "20260525"), 20260525L)
  nstart_mle <- parse_int0(arg0(opts, "nstart-mle", "12"), 12L)
  nstart_gmm <- parse_int0(arg0(opts, "nstart-gmm", "10"), 10L)
  nstart_aug <- parse_int0(arg0(opts, "nstart-aug", as.character(nstart_gmm)), nstart_gmm)
  rho_max <- as.numeric(arg0(opts, "rho-max", "0.98"))
  run_gmm <- parse_bool0(arg0(opts, "run-gmm", "TRUE"), TRUE)
  run_aug <- parse_bool0(arg0(opts, "run-aug", "TRUE"), TRUE)
  aug_n_grid <- parse_int0(arg0(opts, "aug-n-grid", "0"), 0L) # legacy: converted to hmax = T / n_grid when no hmax is supplied
  if (is.na(aug_n_grid) || aug_n_grid <= 0L) aug_n_grid <- NULL
  aug_quad_hmax <- parse_num0(arg0(opts, "aug-quad-hmax", ""), NA_real_)
  if (!is.finite(aug_quad_hmax) || aug_quad_hmax <= 0) aug_quad_hmax <- NULL
  aug_degree <- parse_int0(arg0(opts, "aug-degree", "1"), 1L)
  aug_s <- parse_num_vec0(arg0(opts, "aug-s", "0.4,0.4"), c(0.4, 0.4))
  aug_ridge_rel <- parse_num0(arg0(opts, "aug-ridge-rel", "1e-8"), 1e-8)
  compute_info <- parse_bool0(arg0(opts, "compute-info", "FALSE"), FALSE)
  truth_start <- parse_bool0(arg0(opts, "truth-start", "FALSE"), FALSE)
  overwrite <- parse_bool0(arg0(opts, "overwrite", "FALSE"), FALSE)
  dry_run <- parse_bool0(arg0(opts, "dry-run", "FALSE"), FALSE)
  tmp_build <- parse_bool0(arg0(opts, "tmp-build", "TRUE"), TRUE)

  if (array_id < 1L || array_id > array_count) stop("array-id must be between 1 and array-count")
  raw_file <- file.path(outdir, sprintf("raw_array%04d_of%04d.csv", array_id, array_count))
  fail_file <- file.path(outdir, sprintf("failures_array%04d_of%04d.csv", array_id, array_count))
  manifest_file <- file.path(outdir, sprintf("manifest_array%04d_of%04d.csv", array_id, array_count))
  if (file.exists(raw_file) && !overwrite) { message("Estimator output exists and overwrite=FALSE: ", raw_file); return(invisible(NULL)) }

  load_core(base_dir, tmp_build)
  tasks <- do.call(rbind, lapply(T_grid, function(TT) data.frame(T = TT, rep = seq_len(reps_per_T), stringsAsFactors = FALSE)))
  tasks$task_index <- seq_len(nrow(tasks))
  tasks$array_id <- ((tasks$task_index - 1L) %% array_count) + 1L
  tasks$seed <- vapply(seq_len(nrow(tasks)), function(i) stable_seed(seed, paste0(bivar_spec()$name, "_estimator"), tasks$T[i], tasks$rep[i]), integer(1))
  if (any(duplicated(tasks$seed))) stop("Duplicate seeds in estimator grid")
  my_tasks <- tasks[tasks$array_id == array_id, , drop = FALSE]
  write.csv(my_tasks, manifest_file, row.names = FALSE)
  message("Estimator array ", array_id, "/", array_count, " has ", nrow(my_tasks), " tasks")
  if (dry_run) return(invisible(my_tasks))

  one <- function(ii) {
    row <- my_tasks[ii, ]
    message(sprintf("estimator array %d/%d task %d/%d: T=%s rep=%s", array_id, array_count, ii, nrow(my_tasks), row$T, row$rep))
    tryCatch({
      ans <- run_one_rep(T = as.numeric(row$T), rep_id = as.integer(row$rep), seed = as.integer(row$seed),
                         nstart_mle = nstart_mle, nstart_gmm = nstart_gmm,
                         nstart_aug = nstart_aug, run_gmm = run_gmm, run_aug = run_aug,
                         rho_max = rho_max, compute_info = compute_info, truth_start = truth_start,
                         aug_n_grid = aug_n_grid, aug_quad_hmax = aug_quad_hmax,
                         aug_degree = aug_degree, s_aug = aug_s,
                         aug_ridge_rel = aug_ridge_rel)
      ans$task_index <- row$task_index; ans$array_id <- array_id; ans$seed <- row$seed
      list(ok = TRUE, data = ans, failure = NULL)
    }, error = function(e) {
      list(ok = FALSE, data = NULL,
           failure = data.frame(T = row$T, rep = row$rep, task_index = row$task_index,
                                array_id = array_id, seed = row$seed,
                                error = conditionMessage(e), stringsAsFactors = FALSE))
    })
  }

  n <- nrow(my_tasks)
  if (n == 0L) { write.csv(data.frame(), raw_file, row.names = FALSE); return(invisible(NULL)) }
  workers <- max(1L, min(as.integer(workers), n))
  res <- if (workers > 1L && .Platform$OS.type != "windows") parallel::mclapply(seq_len(n), one, mc.cores = workers, mc.preschedule = FALSE, mc.set.seed = FALSE) else lapply(seq_len(n), one)
  ok <- lapply(res, `[[`, "data"); ok <- ok[!vapply(ok, is.null, logical(1))]
  fail <- lapply(res, `[[`, "failure"); fail <- fail[!vapply(fail, is.null, logical(1))]
  if (length(ok)) write.csv(do.call(rbind, ok), raw_file, row.names = FALSE) else write.csv(data.frame(), raw_file, row.names = FALSE)
  if (length(fail)) write.csv(do.call(rbind, fail), fail_file, row.names = FALSE)
  invisible(NULL)
}

main()

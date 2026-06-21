################################################################################
# bivar_core.R
# Focused bivariate Hawkes simulation utilities for Godambe optimality.
################################################################################

suppressPackageStartupMessages({ library(Rcpp) })
sourceCpp("hawkes_bivar_fast.cpp")

param_names <- function() c("mu1", "mu2", "alpha11", "alpha12", "alpha21", "alpha22", "beta")
param_labels <- function() {
  c(mu1 = expression(mu[1]), mu2 = expression(mu[2]),
    alpha11 = expression(alpha[11]), alpha12 = expression(alpha[12]),
    alpha21 = expression(alpha[21]), alpha22 = expression(alpha[22]), beta = expression(beta))
}

method_order <- function() c("MLE", "GMM_Dtheta", "GMM_aug")
method_labels <- function() {
  c(MLE = "MLE",
    GMM_Dtheta = "just-identified Dtheta GMM",
    GMM_aug = "overidentified bounded-inverse GMM")
}
aug_default_degree <- function() 1L
aug_degree_sanitize <- function(x = aug_default_degree()) {
  x <- as.integer(x)
  if (!is.finite(x) || is.na(x)) x <- aug_default_degree()
  max(1L, min(2L, x))
}
aug_q <- function(degree = aug_default_degree()) 7L * (aug_degree_sanitize(degree) + 1L)
aug_moment_names <- function(degree = aug_default_degree(), q = NULL) {
  if (!is.null(q) && is.finite(q) && q %% length(param_names()) == 0L) {
    degree <- as.integer(q / length(param_names()) - 1L)
  }
  degree <- aug_degree_sanitize(degree)
  blocks <- c("Dtheta", "binv1", "binv2")[seq_len(degree + 1L)]
  paste(rep(blocks, each = length(param_names())),
        rep(param_names(), times = length(blocks)), sep = ":")
}

param_group <- function(p) {
  out <- rep("other", length(p)); names(out) <- p
  out[p %in% c("mu1", "mu2")] <- "baselines"
  out[p %in% c("alpha11", "alpha22")] <- "self-excitation"
  out[p %in% c("alpha12", "alpha21")] <- "cross-excitation"
  out[p %in% "beta"] <- "decay"
  out
}

par_to_vec <- function(mu1, mu2, alpha11, alpha12, alpha21, alpha22, beta) {
  c(mu1 = mu1, mu2 = mu2, alpha11 = alpha11, alpha12 = alpha12,
    alpha21 = alpha21, alpha22 = alpha22, beta = beta)
}

vec_to_alpha <- function(theta) {
  matrix(c(theta[["alpha11"]], theta[["alpha21"]],
           theta[["alpha12"]], theta[["alpha22"]]), nrow = 2, ncol = 2,
         dimnames = list(c("recv1", "recv2"), c("src1", "src2")))
}

spectral_radius_alpha <- function(alpha) {
  ev <- eigen(matrix(alpha, 2, 2), only.values = TRUE)$values
  max(Mod(ev))
}

stationary_rates <- function(mu, alpha) {
  as.numeric(solve(diag(2) - matrix(alpha, 2, 2), as.numeric(mu)))
}

bivar_spec <- function() {
  theta <- par_to_vec(mu1 = 0.22, mu2 = 0.18,
                      alpha11 = 0.34, alpha12 = 0.10,
                      alpha21 = 0.24, alpha22 = 0.30,
                      beta = 1.25)
  alpha <- vec_to_alpha(theta)
  mu <- theta[c("mu1", "mu2")]
  rates <- stationary_rates(mu, alpha)
  list(
    name = "moderate_bivariate_shared_beta",
    A = 3,
    theta = theta,
    mu = mu,
    alpha = alpha,
    beta = theta[["beta"]],
    rho = spectral_radius_alpha(alpha),
    lambda_bar = rates,
    description = "Bivariate normalized truncated-exponential Hawkes with asymmetric cross-excitation and one shared decay"
  )
}

bivar_bounds <- function() {
  list(
    lower = c(mu1 = 0.03, mu2 = 0.03,
              alpha11 = 0.02, alpha12 = 0.02, alpha21 = 0.02, alpha22 = 0.02,
              beta = 0.25),
    upper = c(mu1 = 1.30, mu2 = 1.30,
              alpha11 = 0.75, alpha12 = 0.60, alpha21 = 0.60, alpha22 = 0.75,
              beta = 4.50)
  )
}

burnin_length <- function(A, rho) max(150 * A, 40 * A / (1 - rho))

stable_seed <- function(master_seed, scenario, T, rep_id) {
  m <- 2147483647
  h <- as.numeric(abs(as.integer(master_seed))) %% m
  vals <- c(utf8ToInt(as.character(scenario)), as.integer(round(as.numeric(T) * 1000)), as.integer(rep_id))
  for (z in vals) h <- (h * 16807 + as.numeric(z) + 1) %% m
  seed <- as.integer(h)
  if (seed <= 0L) seed <- 1L
  seed
}

parse_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list(); i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (!startsWith(a, "--")) { i <- i + 1L; next }
    a <- substring(a, 3L)
    if (grepl("=", a, fixed = TRUE)) {
      out[[sub("=.*$", "", a)]] <- sub("^[^=]*=", "", a)
      i <- i + 1L
    } else {
      key <- a
      if (i < length(args) && !startsWith(args[[i + 1L]], "--")) {
        out[[key]] <- args[[i + 1L]]; i <- i + 2L
      } else { out[[key]] <- "TRUE"; i <- i + 1L }
    }
  }
  out
}
script_dir <- function() {
  ff <- commandArgs(FALSE); ff <- ff[grepl("^--file=", ff)]
  if (length(ff) > 0L) return(dirname(normalizePath(sub("^--file=", "", ff[1L]))))
  getwd()
}
arg_value <- function(opts, name, default = NULL, env = character()) {
  if (!is.null(opts[[name]])) return(opts[[name]])
  for (e in env) { v <- Sys.getenv(e, unset = ""); if (nzchar(v)) return(v) }
  default
}
parse_int <- function(x, default) if (is.null(x) || !nzchar(as.character(x))) default else as.integer(x)
parse_num <- function(x, default) if (is.null(x) || !nzchar(as.character(x))) default else as.numeric(x)
parse_bool <- function(x, default = FALSE) {
  if (is.null(x) || !nzchar(as.character(x))) return(default)
  tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y", "on")
}
parse_num_vec <- function(x, default) {
  if (is.null(x) || !nzchar(as.character(x))) return(default)
  z <- suppressWarnings(as.numeric(strsplit(as.character(x), ",")[[1L]]))
  z <- z[is.finite(z)]
  if (!length(z)) default else z
}

aug_default_hmax <- function(T, A = NULL) {
  # Event/expiration-adaptive Gauss-Legendre integration is the default.  The
  # hmax cap only splits long smooth intervals; discontinuities are always
  # split exactly at event and support-expiration times.
  if (!is.null(A) && is.finite(A) && A > 0) return(min(0.50, max(0.10, A / 6)))
  0.50
}
aug_hmax_from_legacy_grid <- function(T, n_grid = NULL, A = NULL) {
  if (!is.null(n_grid) && is.finite(n_grid) && n_grid > 0) return(as.numeric(T) / as.numeric(n_grid))
  aug_default_hmax(T, A)
}
solve_spd_unregularized <- function(M) {
  M <- symmetrize(as.matrix(M))
  ch <- try(chol(M), silent = TRUE)
  if (!inherits(ch, "try-error")) return(chol2inv(ch))
  solve(M)
}
regularized_weight_from_omega <- function(Omega, ridge_rel = 1e-8, tol = 1e-10,
                                          cond_max = 1e12, row_scale = TRUE) {
  Omega <- symmetrize(as.matrix(Omega))
  q <- nrow(Omega)
  if (q != ncol(Omega)) stop("Omega must be square")
  d <- rep(1, q)
  if (isTRUE(row_scale)) {
    d0 <- sqrt(pmax(diag(Omega), .Machine$double.eps))
    d0[!is.finite(d0) | d0 <= 0] <- 1
    d <- 1 / d0
  }
  S <- diag(d, q)
  Omega_s <- symmetrize(S %*% Omega %*% S)
  ev <- try(eigen(Omega_s, symmetric = TRUE, only.values = TRUE)$values, silent = TRUE)
  min_eig <- if (inherits(ev, "try-error") || !length(ev)) NA_real_ else min(ev)
  max_eig <- if (inherits(ev, "try-error") || !length(ev)) NA_real_ else max(ev)
  cond <- if (is.finite(min_eig) && min_eig > 0 && is.finite(max_eig)) max_eig / min_eig else Inf
  activated <- !is.finite(min_eig) || min_eig <= tol || !is.finite(cond) || cond > cond_max
  ridge_abs_scaled <- 0
  Omega_used_s <- Omega_s
  if (activated) {
    scale_s <- mean(diag(Omega_s))
    if (!is.finite(scale_s) || scale_s <= 0) scale_s <- 1
    ridge_abs_scaled <- ridge_rel * scale_s
    Omega_used_s <- Omega_s + diag(ridge_abs_scaled, q)
  }
  W_s <- try(solve_spd_unregularized(Omega_used_s), silent = TRUE)
  if (inherits(W_s, "try-error") || any(!is.finite(W_s))) {
    activated <- TRUE
    scale_s <- mean(diag(Omega_s))
    if (!is.finite(scale_s) || scale_s <= 0) scale_s <- 1
    ridge_abs_scaled <- max(ridge_abs_scaled, ridge_rel * scale_s, tol * scale_s)
    Omega_used_s <- Omega_s + diag(ridge_abs_scaled, q)
    W_s <- solve_spd_unregularized(Omega_used_s)
  }
  W <- symmetrize(S %*% W_s %*% S)
  list(W = W, Omega = Omega, Omega_scaled = Omega_s, Omega_used_scaled = Omega_used_s,
       ridge_active = activated, ridge_abs = ridge_abs_scaled, min_eig = min_eig,
       max_eig = max_eig, cond = cond, row_scale = d)
}


valid_theta <- function(theta, lower, upper, rho_max = 0.98) {
  theta <- as.numeric(theta); names(theta) <- names(lower)
  if (length(theta) != 7L || any(!is.finite(theta))) return(FALSE)
  if (any(theta < lower) || any(theta > upper)) return(FALSE)
  if (theta[["mu1"]] <= 0 || theta[["mu2"]] <= 0 || theta[["beta"]] <= 0) return(FALSE)
  if (any(theta[c("alpha11", "alpha12", "alpha21", "alpha22")] <= 0)) return(FALSE)
  spectral_radius_alpha(vec_to_alpha(theta)) < rho_max
}
clamp <- function(x, lower, upper) { x <- pmin(pmax(x, lower), upper); names(x) <- names(lower); x }

random_theta <- function(lower, upper, rho_max = 0.98, max_tries = 10000L) {
  for (i in seq_len(max_tries)) {
    th <- runif(7, lower, upper); names(th) <- names(lower)
    if (valid_theta(th, lower, upper, rho_max)) return(th)
  }
  stop("Could not generate stable random bivariate start")
}

moment_start <- function(ev, T, bounds, rho_max = 0.98) {
  rate_hat <- c(sum(ev$type == 1L & ev$time >= 0 & ev$time <= T),
                sum(ev$type == 2L & ev$time >= 0 & ev$time <= T)) / T
  rate_hat <- pmax(rate_hat, 1e-4)
  candidates_alpha <- list(
    matrix(c(0.30, 0.20, 0.08, 0.28), 2, 2),
    matrix(c(0.34, 0.24, 0.10, 0.30), 2, 2),
    matrix(c(0.25, 0.15, 0.15, 0.25), 2, 2),
    matrix(c(0.42, 0.16, 0.06, 0.26), 2, 2),
    matrix(c(0.25, 0.28, 0.12, 0.33), 2, 2)
  )
  betas <- c(0.80, 1.25, 1.80)
  rows <- list()
  for (A0 in candidates_alpha) {
    if (spectral_radius_alpha(A0) >= rho_max) next
    mu0 <- as.numeric((diag(2) - A0) %*% rate_hat)
    mu0 <- pmax(mu0, bounds$lower[c("mu1", "mu2")])
    for (b in betas) {
      th <- c(mu1 = mu0[1], mu2 = mu0[2],
              alpha11 = A0[1, 1], alpha12 = A0[1, 2],
              alpha21 = A0[2, 1], alpha22 = A0[2, 2], beta = b)
      th <- clamp(th, bounds$lower, bounds$upper)
      if (valid_theta(th, bounds$lower, bounds$upper, rho_max)) rows[[length(rows) + 1L]] <- th
    }
  }
  if (!length(rows)) return(matrix(numeric(0), 0, 7, dimnames = list(NULL, names(bounds$lower))))
  do.call(rbind, rows)
}

make_starts <- function(lower, upper, nstart, theta_hint = NULL, rho_max = 0.98,
                        ev = NULL, T = NULL, include_truth = FALSE) {
  starts <- list()
  add_start <- function(th) {
    th <- as.numeric(th); names(th) <- names(lower)
    th <- clamp(th, lower, upper)
    if (valid_theta(th, lower, upper, rho_max)) starts[[length(starts) + 1L]] <<- th
  }
  if (!is.null(theta_hint)) {
    add_start(theta_hint)
    for (sd in c(0.05, 0.12, 0.25)) add_start(theta_hint * exp(rnorm(7, sd = sd)))
  }
  if (!is.null(ev) && !is.null(T)) {
    ms <- moment_start(ev, T, list(lower = lower, upper = upper), rho_max = rho_max)
    if (nrow(ms)) for (i in seq_len(nrow(ms))) add_start(ms[i, ])
  }
  while (length(starts) < nstart) starts[[length(starts) + 1L]] <- random_theta(lower, upper, rho_max)
  mat <- do.call(rbind, starts)
  keep <- !duplicated(round(mat, 10))
  starts <- lapply(seq_len(nrow(mat))[keep], function(i) { z <- mat[i, ]; names(z) <- names(lower); z })
  while (length(starts) < nstart) starts[[length(starts) + 1L]] <- random_theta(lower, upper, rho_max)
  starts[seq_len(nstart)]
}

simulate_bivar <- function(T, spec = bivar_spec(), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  simulate_hawkes_bivar_cluster_cpp(T, spec$A, as.numeric(spec$mu), spec$alpha,
                                    spec$beta, burnin_length(spec$A, spec$rho))
}

loglik_score <- function(theta, ev, T, A, rho_max = 0.98) {
  hawkes_bivar_loglik_score_cpp(as.numeric(theta), ev$time, as.integer(ev$type), T, A, rho_max)
}

finite_diff_info <- function(theta, ev, T, A, rho_max = 0.98) {
  p <- 7L
  H <- matrix(NA_real_, p, p)
  score_at <- function(th) as.numeric(loglik_score(th, ev, T, A, rho_max)$score)
  for (k in seq_len(p)) {
    h <- 1e-5 * max(1, abs(theta[k]))
    thp <- theta; thm <- theta
    thp[k] <- thp[k] + h
    thm[k] <- thm[k] - h
    sp <- score_at(thp); sm <- score_at(thm)
    H[, k] <- (sp - sm) / (2 * h)
  }
  info <- -(H + t(H)) / 2
  rownames(info) <- colnames(info) <- param_names()
  info
}

fit_mle <- function(ev, T, A, bounds = bivar_bounds(), nstart = 12,
                    theta_hint = NULL, rho_max = 0.98, compute_info = FALSE,
                    maxit = 1200) {
  lower <- bounds$lower; upper <- bounds$upper
  starts <- make_starts(lower, upper, nstart, theta_hint = theta_hint, rho_max = rho_max, ev = ev, T = T)
  best <- NULL
  for (st in starts) {
    cache_theta <- NULL; cache_res <- NULL
    eval <- function(th) {
      if (is.null(cache_theta) || any(th != cache_theta)) {
        cache_theta <<- th; cache_res <<- loglik_score(th, ev, T, A, rho_max)
      }
      cache_res
    }
    fn <- function(th) {
      r <- eval(th)
      if (!isTRUE(r$valid) || !is.finite(r$loglik)) return(1e100)
      -as.numeric(r$loglik)
    }
    gr <- function(th) {
      r <- eval(th)
      if (!isTRUE(r$valid) || !all(is.finite(r$score))) return(rep(0, 7))
      -as.numeric(r$score)
    }
    fit <- try(optim(st, fn, gr, method = "L-BFGS-B", lower = lower, upper = upper,
                     control = list(maxit = maxit, factr = 1e7)), silent = TRUE)
    if (!inherits(fit, "try-error") && is.finite(fit$value)) {
      if (is.null(best) || fit$value < best$value) best <- fit
    }
  }
  if (is.null(best)) stop("All MLE starts failed")
  names(best$par) <- names(lower)
  info <- vcov <- matrix(NA_real_, 7, 7, dimnames = list(param_names(), param_names()))
  pd <- FALSE; min_info_eig <- NA_real_
  if (compute_info) {
    it <- try(finite_diff_info(best$par, ev, T, A, rho_max), silent = TRUE)
    if (!inherits(it, "try-error") && all(is.finite(it))) {
      info <- it
      evs <- eigen((info + t(info)) / 2, symmetric = TRUE, only.values = TRUE)$values
      min_info_eig <- min(evs)
      pd <- is.finite(min_info_eig) && min_info_eig > 1e-8
      if (pd) vcov <- solve(info)
    }
  }
  list(par = best$par, value = best$value, converged = best$convergence == 0,
       info = info, vcov = vcov, info_pd = pd, min_info_eig = min_info_eig,
       optim_message = if (!is.null(best$message)) best$message else "")
}

fit_gmm_dtheta <- function(ev, T, A, bounds = bivar_bounds(), nstart = 10,
                           theta_hint = NULL, rho_max = 0.98, maxit = 1000) {
  lower <- bounds$lower
  upper <- bounds$upper

  starts <- make_starts(
    lower, upper, nstart,
    theta_hint = theta_hint,
    rho_max = rho_max,
    ev = ev,
    T = T
  )

  eps <- 1e-8

  to_z <- function(th) {
    th <- as.numeric(th)
    names(th) <- names(lower)
    u <- (th - lower) / (upper - lower)
    u <- pmin(pmax(u, eps), 1 - eps)
    qlogis(u)
  }

  from_z <- function(z) {
    th <- lower + (upper - lower) * plogis(z)
    names(th) <- names(lower)
    th
  }

  mom <- function(th) {
    hawkes_bivar_ls_moment_exact_cpp(
      as.numeric(th),
      ev$time,
      as.integer(ev$type),
      T,
      A,
      rho_max
    )
  }

  obj_theta <- function(th) {
    th <- as.numeric(th)
    names(th) <- names(lower)

    if (!valid_theta(th, lower, upper, rho_max)) return(1e100)

    m <- mom(th)
    if (any(!is.finite(m))) return(1e100)

    # Same minimizer as sum(m^2), but numerically better because m is O(T^{-1/2}).
    as.numeric(T) * sum(m * m)
  }

  obj_z <- function(z) {
    obj_theta(from_z(z))
  }

  best <- NULL

  consider <- function(par, value, convergence = NA_integer_, message = "") {
    if (!is.finite(value)) return()
    th <- from_z(par)
    val <- obj_theta(th)
    if (!is.finite(val)) return()
    fit <- list(par = th, value = val, convergence = convergence, message = message)
    if (is.null(best) || fit$value < best$value) best <<- fit
  }

  for (st in starts) {
    z0 <- to_z(st)

    # Always include the starting point itself.
    consider(z0, obj_z(z0), convergence = 999L, message = "start")

    # BFGS on transformed parameters avoids L-BFGS-B falsely accepting the MLE start.
    fit1 <- try(
      optim(
        z0,
        obj_z,
        method = "BFGS",
        control = list(maxit = maxit, reltol = 1e-12)
      ),
      silent = TRUE
    )

    if (!inherits(fit1, "try-error") && is.finite(fit1$value)) {
      consider(
        fit1$par,
        fit1$value,
        convergence = fit1$convergence,
        message = if (!is.null(fit1$message)) fit1$message else "BFGS"
      )
    }

    # A short Nelder-Mead fallback helps when finite-difference BFGS is sticky.
    fit2 <- try(
      optim(
        z0,
        obj_z,
        method = "Nelder-Mead",
        control = list(maxit = max(300L, floor(maxit / 3L)), reltol = 1e-10)
      ),
      silent = TRUE
    )

    if (!inherits(fit2, "try-error") && is.finite(fit2$value)) {
      consider(
        fit2$par,
        fit2$value,
        convergence = fit2$convergence,
        message = if (!is.null(fit2$message)) fit2$message else "Nelder-Mead"
      )
    }
  }

  if (is.null(best)) stop("All GMM starts failed")

  names(best$par) <- names(lower)

  list(
    par = best$par,
    value = best$value,
    converged = best$convergence == 0,
    moment_norm = sqrt(best$value / as.numeric(T)),
    optim_message = if (!is.null(best$message)) best$message else ""
  )
}




fit_gmm_aug <- function(ev, T, A, bounds = bivar_bounds(), nstart = 10,
                        theta_hint = NULL, rho_max = 0.98, maxit = 1000,
                        n_grid = NULL, aug_quad_hmax = NULL,
                        s_aug = c(0.4, 0.4), aug_degree = aug_default_degree(),
                        ridge_rel = 1e-8, weight_tol = 1e-10,
                        weight_cond_max = 1e12) {
  lower <- bounds$lower
  upper <- bounds$upper
  s_aug <- as.numeric(s_aug)
  if (length(s_aug) == 1L) s_aug <- rep(s_aug, 2L)
  if (length(s_aug) != 2L || any(!is.finite(s_aug)) || any(s_aug <= 0)) {
    stop("s_aug must contain one or two positive finite stabilizers")
  }
  aug_degree <- aug_degree_sanitize(aug_degree)
  q_aug <- aug_q(aug_degree)
  if (is.null(aug_quad_hmax) || !is.finite(aug_quad_hmax) || aug_quad_hmax <= 0) {
    aug_quad_hmax <- aug_hmax_from_legacy_grid(T, n_grid = n_grid, A = A)
  }
  aug_quad_hmax <- as.numeric(aug_quad_hmax)

  starts <- make_starts(
    lower, upper, nstart,
    theta_hint = theta_hint,
    rho_max = rho_max,
    ev = ev,
    T = T
  )

  eps <- 1e-8

  to_z <- function(th) {
    th <- as.numeric(th)
    names(th) <- names(lower)
    u <- (th - lower) / (upper - lower)
    u <- pmin(pmax(u, eps), 1 - eps)
    qlogis(u)
  }

  from_z <- function(z) {
    th <- lower + (upper - lower) * plogis(z)
    names(th) <- names(lower)
    th
  }

  mom <- function(th) {
    hawkes_bivar_aug_moment_adaptive_cpp(
      as.numeric(th), ev$time, as.integer(ev$type), T, A,
      rho_max, s_aug[1], s_aug[2], aug_degree, aug_quad_hmax
    )
  }

  obj_theta <- function(th, W = NULL) {
    th <- as.numeric(th)
    names(th) <- names(lower)
    if (!valid_theta(th, lower, upper, rho_max)) return(1e100)
    m <- mom(th)
    if (length(m) != q_aug || any(!is.finite(m))) return(1e100)
    if (is.null(W)) return(as.numeric(T) * sum(m * m))
    val <- as.numeric(T) * drop(crossprod(m, W %*% m))
    if (!is.finite(val)) 1e100 else val
  }

  optimize_stage <- function(stage_starts, W = NULL, label = "stage", maxit_stage = maxit) {
    best <- NULL
    consider <- function(z, value, convergence = NA_integer_, message = "") {
      if (!is.finite(value)) return()
      th <- from_z(z)
      val <- obj_theta(th, W)
      if (!is.finite(val)) return()
      m <- mom(th)
      fit <- list(par = th, value = val, convergence = convergence,
                  message = message, moment = m)
      if (is.null(best) || fit$value < best$value) best <<- fit
    }
    for (st in stage_starts) {
      z0 <- to_z(st)
      obj_z <- function(z) obj_theta(from_z(z), W)
      consider(z0, obj_z(z0), convergence = 999L, message = paste(label, "start"))
      fit1 <- try(
        optim(z0, obj_z, method = "BFGS",
              control = list(maxit = maxit_stage, reltol = 1e-12)),
        silent = TRUE
      )
      if (!inherits(fit1, "try-error") && is.finite(fit1$value)) {
        consider(fit1$par, fit1$value, convergence = fit1$convergence,
                 message = if (!is.null(fit1$message)) fit1$message else paste(label, "BFGS"))
      }
      fit2 <- try(
        optim(z0, obj_z, method = "Nelder-Mead",
              control = list(maxit = max(300L, floor(maxit_stage / 3L)), reltol = 1e-10)),
        silent = TRUE
      )
      if (!inherits(fit2, "try-error") && is.finite(fit2$value)) {
        consider(fit2$par, fit2$value, convergence = fit2$convergence,
                 message = if (!is.null(fit2$message)) fit2$message else paste(label, "Nelder-Mead"))
      }
    }
    best
  }

  first <- optimize_stage(starts, W = NULL, label = "aug first", maxit_stage = maxit)
  if (is.null(first)) stop("All augmented GMM first-step starts failed")

  om <- hawkes_bivar_aug_omega_adaptive_cpp(
    as.numeric(first$par), ev$time, as.integer(ev$type), T, A,
    rho_max, s_aug[1], s_aug[2], aug_degree, aug_quad_hmax
  )
  if (!isTRUE(om$valid) || as.integer(om$n_eval) <= 0L) stop("Augmented GMM Omega estimate failed")
  reg <- regularized_weight_from_omega(om$Omega_aug, ridge_rel = ridge_rel,
                                       tol = weight_tol, cond_max = weight_cond_max,
                                       row_scale = TRUE)

  second_starts <- c(list(first$par), starts)
  if (!is.null(theta_hint)) second_starts <- c(list(theta_hint), second_starts)
  smat <- do.call(rbind, lapply(second_starts, function(x) { x <- as.numeric(x); names(x) <- names(lower); clamp(x, lower, upper) }))
  smat <- smat[!duplicated(round(smat, 10)), , drop = FALSE]
  second_starts <- lapply(seq_len(nrow(smat)), function(i) { z <- smat[i, ]; names(z) <- names(lower); z })

  second <- optimize_stage(second_starts, W = reg$W, label = "aug second", maxit_stage = maxit)
  if (is.null(second)) stop("All augmented GMM second-step starts failed")
  names(second$par) <- names(lower)
  m_final <- mom(second$par)

  list(
    par = second$par,
    value = second$value,
    converged = second$convergence == 0,
    moment_norm = sqrt(sum(m_final * m_final)),
    weighted_moment_norm = sqrt(max(0, drop(crossprod(m_final, reg$W %*% m_final)))),
    stage1_value = first$value,
    stage1_moment_norm = sqrt(sum(first$moment * first$moment)),
    stage2_value = second$value,
    ridge_active = isTRUE(reg$ridge_active),
    ridge_abs = reg$ridge_abs,
    omega_min_eig = reg$min_eig,
    omega_cond = reg$cond,
    omega_n_eval = as.integer(om$n_eval),
    omega_n_intervals = as.integer(om$n_intervals),
    n_grid = as.integer(om$n_eval),
    aug_quad_hmax = aug_quad_hmax,
    aug_degree = aug_degree,
    q_aug = q_aug,
    s_aug1 = s_aug[1],
    s_aug2 = s_aug[2],
    optim_message = if (!is.null(second$message)) second$message else ""
  )
}



one_estimate_to_long <- function(theta_hat, theta_true, method, T, rep_id,
                                 converged = TRUE, vcov = NULL, boundary_hit = FALSE,
                                 info_pd = NA, criterion = NA, min_info_eig = NA,
                                 lower = NULL, upper = NULL, boundary_tol = 1e-5,
                                 ridge_active = NA, ridge_abs = NA,
                                 weighted_criterion = NA,
                                 stage1_moment_norm = NA,
                                 weighted_moment_norm = NA,
                                 omega_min_eig = NA,
                                 omega_cond = NA,
                                 omega_n_eval = NA,
                                 omega_n_intervals = NA,
                                 aug_quad_hmax = NA,
                                 aug_degree = NA,
                                 q_aug = NA,
                                 n_grid = NA) {
  se <- rep(NA_real_, length(theta_hat))
  if (!is.null(vcov) && all(dim(vcov) == c(length(theta_hat), length(theta_hat))) && all(is.finite(diag(vcov)))) {
    se <- sqrt(pmax(diag(vcov), 0))
  }
  if (is.null(lower)) lower <- rep(NA_real_, length(theta_hat))
  if (is.null(upper)) upper <- rep(NA_real_, length(theta_hat))
  at_lower <- is.finite(lower) & abs(theta_hat - lower) <= boundary_tol
  at_upper <- is.finite(upper) & abs(theta_hat - upper) <= boundary_tol
  data.frame(
    scenario = bivar_spec()$name,
    T = T,
    rep = rep_id,
    method = method,
    param = names(theta_true),
    estimate = as.numeric(theta_hat),
    true = as.numeric(theta_true),
    error = as.numeric(theta_hat - theta_true),
    se = se,
    cover95 = ifelse(is.finite(se), abs(theta_hat - theta_true) <= 1.96 * se, NA),
    converged = converged,
    boundary_hit = boundary_hit,
    at_lower = as.logical(at_lower),
    at_upper = as.logical(at_upper),
    lower = as.numeric(lower),
    upper = as.numeric(upper),
    info_pd = info_pd,
    min_info_eig = min_info_eig,
    criterion = criterion,
    ridge_active = ridge_active,
    ridge_abs = ridge_abs,
    weighted_criterion = weighted_criterion,
    stage1_moment_norm = stage1_moment_norm,
    weighted_moment_norm = weighted_moment_norm,
    omega_min_eig = omega_min_eig,
    omega_cond = omega_cond,
    omega_n_eval = omega_n_eval,
    omega_n_intervals = omega_n_intervals,
    aug_quad_hmax = aug_quad_hmax,
    aug_degree = aug_degree,
    q_aug = q_aug,
    n_grid = n_grid,
    stringsAsFactors = FALSE
  )
}

run_one_rep <- function(T, rep_id, spec = bivar_spec(), bounds = bivar_bounds(),
                        nstart_mle = 12, nstart_gmm = 10, nstart_aug = nstart_gmm,
                        run_gmm = TRUE, run_aug = TRUE,
                        rho_max = 0.98, seed = NULL, compute_info = FALSE,
                        truth_start = FALSE, aug_n_grid = NULL,
                        aug_quad_hmax = NULL, aug_degree = aug_default_degree(),
                        s_aug = c(0.4, 0.4), aug_ridge_rel = 1e-8) {
  ev <- simulate_bivar(T, spec, seed)
  theta_true <- spec$theta
  lower <- bounds$lower; upper <- bounds$upper; tol <- 1e-5
  hint <- if (truth_start) theta_true else NULL
  mle <- fit_mle(ev, T, spec$A, bounds, nstart = nstart_mle, theta_hint = hint,
                 rho_max = rho_max, compute_info = compute_info)
  mle_boundary <- any(abs(mle$par - lower) < tol | abs(mle$par - upper) < tol)
  rows <- list(one_estimate_to_long(mle$par, theta_true, "MLE", T, rep_id,
                                    mle$converged, mle$vcov, mle_boundary, mle$info_pd,
                                    mle$value, mle$min_info_eig, lower, upper))
  if (isTRUE(run_gmm)) {
    gmm <- fit_gmm_dtheta(ev, T, spec$A, bounds, nstart = nstart_gmm,
                          theta_hint = mle$par, rho_max = rho_max)
    gmm_boundary <- any(abs(gmm$par - lower) < tol | abs(gmm$par - upper) < tol)
    rows[[length(rows) + 1L]] <- one_estimate_to_long(gmm$par, theta_true, "GMM_Dtheta", T, rep_id,
                                                       gmm$converged, NULL, gmm_boundary, NA,
                                                       gmm$moment_norm, NA, lower, upper)
  }
  if (isTRUE(run_aug)) {
    aug <- fit_gmm_aug(ev, T, spec$A, bounds, nstart = nstart_aug,
                       theta_hint = mle$par, rho_max = rho_max,
                       n_grid = aug_n_grid, aug_quad_hmax = aug_quad_hmax,
                       aug_degree = aug_degree, s_aug = s_aug,
                       ridge_rel = aug_ridge_rel)
    aug_boundary <- any(abs(aug$par - lower) < tol | abs(aug$par - upper) < tol)
    rows[[length(rows) + 1L]] <- one_estimate_to_long(
      aug$par, theta_true, "GMM_aug", T, rep_id,
      aug$converged, NULL, aug_boundary, NA,
      aug$moment_norm, NA, lower, upper,
      ridge_active = aug$ridge_active,
      ridge_abs = aug$ridge_abs,
      weighted_criterion = aug$weighted_moment_norm,
      stage1_moment_norm = aug$stage1_moment_norm,
      weighted_moment_norm = aug$weighted_moment_norm,
      omega_min_eig = aug$omega_min_eig,
      omega_cond = aug$omega_cond,
      omega_n_eval = aug$omega_n_eval,
      omega_n_intervals = aug$omega_n_intervals,
      aug_quad_hmax = aug$aug_quad_hmax,
      aug_degree = aug$aug_degree,
      q_aug = aug$q_aug,
      n_grid = aug$n_grid
    )
  }
  ans <- do.call(rbind, rows)
  ans$n_events_total <- sum(ev$time >= 0 & ev$time <= T)
  ans$n_events_1 <- sum(ev$type == 1L & ev$time >= 0 & ev$time <= T)
  ans$n_events_2 <- sum(ev$type == 2L & ev$time >= 0 & ev$time <= T)
  rownames(ans) <- NULL
  ans
}


symmetrize <- function(M) (M + t(M)) / 2
safe_solve <- function(M, tol = 1e-10) {
  M <- symmetrize(as.matrix(M))
  ev <- eigen(M, symmetric = TRUE)
  vals <- pmax(ev$values, tol)
  ev$vectors %*% diag(1 / vals, length(vals)) %*% t(ev$vectors)
}
sym_sqrt <- function(M, inverse = FALSE, tol = 1e-10) {
  M <- symmetrize(as.matrix(M))
  ev <- eigen(M, symmetric = TRUE)
  vals <- pmax(ev$values, tol)
  vals <- if (inverse) 1 / sqrt(vals) else sqrt(vals)
  ev$vectors %*% diag(vals, length(vals)) %*% t(ev$vectors)
}
relative_eigs <- function(V, I, tol = 1e-10) {
  S <- sym_sqrt(I, inverse = FALSE, tol = tol)
  sort(eigen(symmetrize(S %*% V %*% S), symmetric = TRUE, only.values = TRUE)$values)
}

matrix_to_long <- function(M, replicate, matrix_name) {
  M <- as.matrix(M)
  data.frame(scenario = bivar_spec()$name, replicate = replicate, matrix = matrix_name,
             row = rep(seq_len(nrow(M)), times = ncol(M)),
             col = rep(seq_len(ncol(M)), each = nrow(M)),
             value = as.numeric(M), stringsAsFactors = FALSE)
}

godambe_from_matrices <- function(I, A_ls, Omega_ls, A_aug = NULL, Omega_aug = NULL, tol = 1e-10) {
  I <- symmetrize(I)
  A_ls <- symmetrize(A_ls)
  Omega_ls <- symmetrize(Omega_ls)
  I_inv <- safe_solve(I, tol)
  A_inv <- safe_solve(A_ls, tol)
  V_ls <- symmetrize(A_inv %*% Omega_ls %*% A_inv)
  G_ls <- symmetrize(t(A_ls) %*% safe_solve(Omega_ls, tol) %*% A_ls)
  eig_cov <- relative_eigs(V_ls, I, tol)
  S_inv <- sym_sqrt(I, inverse = TRUE, tol = tol)
  eig_info <- sort(eigen(symmetrize(S_inv %*% G_ls %*% S_inv), symmetric = TRUE, only.values = TRUE)$values)
  out <- list(I = I, A_ls = A_ls, Omega_ls = Omega_ls, I_inv = I_inv, V_ls = V_ls, G_ls = G_ls,
              eig_cov_ls_to_mle = eig_cov, eig_info_ls_to_mle = eig_info,
              trace_ratio_ls = sum(diag(V_ls)) / sum(diag(I_inv)),
              logdet_ratio_ls = as.numeric(determinant(V_ls, logarithm = TRUE)$modulus - determinant(I_inv, logarithm = TRUE)$modulus))
  if (!is.null(A_aug) && !is.null(Omega_aug) && length(A_aug) && length(Omega_aug)) {
    A_aug <- as.matrix(A_aug)
    Omega_aug <- symmetrize(Omega_aug)
    G_aug <- symmetrize(t(A_aug) %*% safe_solve(Omega_aug, tol) %*% A_aug)
    V_aug <- safe_solve(G_aug, tol)
    eig_cov_aug <- relative_eigs(V_aug, I, tol)
    eig_info_aug <- sort(eigen(symmetrize(S_inv %*% G_aug %*% S_inv), symmetric = TRUE, only.values = TRUE)$values)
    out$A_aug <- A_aug
    out$Omega_aug <- Omega_aug
    out$V_aug <- V_aug
    out$G_aug <- G_aug
    out$eig_cov_aug_to_mle <- eig_cov_aug
    out$eig_info_aug_to_mle <- eig_info_aug
    out$trace_ratio_aug <- sum(diag(V_aug)) / sum(diag(I_inv))
    out$logdet_ratio_aug <- as.numeric(determinant(V_aug, logarithm = TRUE)$modulus - determinant(I_inv, logarithm = TRUE)$modulus)
  }
  out
}

one_population_task <- function(replicate, T_pop = 80000, n_eval = 80000, seed = 1,
                                spec = bivar_spec(), rho_max = 0.98,
                                s_aug = c(0.4, 0.4),
                                aug_degree = aug_default_degree()) {
  set.seed(seed)
  theta_true <- spec$theta
  ev <- simulate_bivar(T_pop, spec, seed = NULL)
  t_eval <- sort(runif(n_eval, 0, T_pop))
  s_aug <- as.numeric(s_aug)
  if (length(s_aug) == 1L) s_aug <- rep(s_aug, 2L)
  if (length(s_aug) != 2L || any(!is.finite(s_aug)) || any(s_aug <= 0)) {
    stop("s_aug must contain one or two positive finite stabilizers")
  }
  aug_degree <- aug_degree_sanitize(aug_degree)
  mats <- hawkes_bivar_time_average_matrices_cpp(as.numeric(theta_true), ev$time, as.integer(ev$type), t_eval, spec$A,
                                                 rho_max, s_aug[1], s_aug[2], aug_degree)
  if (!isTRUE(mats$valid) || as.integer(mats$n_eval) <= 0L) stop("Matrix time average failed")
  gd <- godambe_from_matrices(mats$I, mats$A_ls, mats$Omega_ls, mats$A_aug, mats$Omega_aug)
  summary <- data.frame(
    scenario = spec$name, replicate = replicate, T_pop = T_pop, n_eval = as.integer(mats$n_eval),
    aug_degree = aug_degree, q_aug = aug_q(aug_degree),
    theta_mu1 = spec$theta[["mu1"]], theta_mu2 = spec$theta[["mu2"]],
    theta_alpha11 = spec$theta[["alpha11"]], theta_alpha12 = spec$theta[["alpha12"]],
    theta_alpha21 = spec$theta[["alpha21"]], theta_alpha22 = spec$theta[["alpha22"]],
    theta_beta = spec$theta[["beta"]],
    rho = spec$rho, lambda_bar1 = spec$lambda_bar[1], lambda_bar2 = spec$lambda_bar[2],
    n_events_total = sum(ev$time >= 0 & ev$time <= T_pop),
    n_events_1 = sum(ev$type == 1L & ev$time >= 0 & ev$time <= T_pop),
    n_events_2 = sum(ev$type == 2L & ev$time >= 0 & ev$time <= T_pop),
    trace_ratio_ls_to_mle = gd$trace_ratio_ls,
    logdet_ratio_ls_to_mle = gd$logdet_ratio_ls,
    eig_min_ls_to_mle = min(gd$eig_cov_ls_to_mle),
    eig_median_ls_to_mle = median(gd$eig_cov_ls_to_mle),
    eig_max_ls_to_mle = max(gd$eig_cov_ls_to_mle),
    trace_ratio_aug_to_mle = if (!is.null(gd$trace_ratio_aug)) gd$trace_ratio_aug else NA_real_,
    logdet_ratio_aug_to_mle = if (!is.null(gd$logdet_ratio_aug)) gd$logdet_ratio_aug else NA_real_,
    eig_min_aug_to_mle = if (!is.null(gd$eig_cov_aug_to_mle)) min(gd$eig_cov_aug_to_mle) else NA_real_,
    eig_median_aug_to_mle = if (!is.null(gd$eig_cov_aug_to_mle)) median(gd$eig_cov_aug_to_mle) else NA_real_,
    eig_max_aug_to_mle = if (!is.null(gd$eig_cov_aug_to_mle)) max(gd$eig_cov_aug_to_mle) else NA_real_,
    stringsAsFactors = FALSE
  )
  matrices <- rbind(matrix_to_long(mats$I, replicate, "I"),
                    matrix_to_long(mats$A_ls, replicate, "A_ls"),
                    matrix_to_long(mats$Omega_ls, replicate, "Omega_ls"),
                    matrix_to_long(mats$A_aug, replicate, "A_aug"),
                    matrix_to_long(mats$Omega_aug, replicate, "Omega_aug"))
  list(summary = summary, matrices = matrices)
}

matrix_from_long <- function(x, matrix_name) {
  xx <- x[x$matrix == matrix_name, , drop = FALSE]
  if (!nrow(xx)) return(NULL)
  nr <- max(xx$row)
  nc <- max(xx$col)
  M <- matrix(0, nr, nc)
  for (i in seq_len(nrow(xx))) M[xx$row[i], xx$col[i]] <- xx$value[i]
  if (nr == 7L) rownames(M) <- param_names() else if (nr %% 7L == 0L) rownames(M) <- aug_moment_names(q = nr)
  if (nc == 7L) colnames(M) <- param_names() else if (nc %% 7L == 0L) colnames(M) <- aug_moment_names(q = nc)
  M
}

build_population_list <- function(mat_long) {
  reps <- sort(unique(mat_long$replicate))
  mats_by_rep <- lapply(reps, function(r) {
    xx <- mat_long[mat_long$replicate == r, , drop = FALSE]
    list(I = matrix_from_long(xx, "I"),
         A_ls = matrix_from_long(xx, "A_ls"),
         Omega_ls = matrix_from_long(xx, "Omega_ls"),
         A_aug = matrix_from_long(xx, "A_aug"),
         Omega_aug = matrix_from_long(xx, "Omega_aug"))
  })
  I_bar <- Reduce(`+`, lapply(mats_by_rep, `[[`, "I")) / length(mats_by_rep)
  A_bar <- Reduce(`+`, lapply(mats_by_rep, `[[`, "A_ls")) / length(mats_by_rep)
  O_bar <- Reduce(`+`, lapply(mats_by_rep, `[[`, "Omega_ls")) / length(mats_by_rep)
  has_aug <- all(vapply(mats_by_rep, function(z) !is.null(z$A_aug) && !is.null(z$Omega_aug), logical(1)))
  if (has_aug) {
    A_aug_bar <- Reduce(`+`, lapply(mats_by_rep, `[[`, "A_aug")) / length(mats_by_rep)
    O_aug_bar <- Reduce(`+`, lapply(mats_by_rep, `[[`, "Omega_aug")) / length(mats_by_rep)
    gd <- godambe_from_matrices(I_bar, A_bar, O_bar, A_aug_bar, O_aug_bar)
  } else {
    gd <- godambe_from_matrices(I_bar, A_bar, O_bar)
  }
  gd
}

population_vcov <- function(pop, method) {
  if (identical(method, "MLE")) return(pop$I_inv)
  if (identical(method, "GMM_Dtheta")) return(pop$V_ls)
  if (identical(method, "GMM_aug") && !is.null(pop$V_aug)) return(pop$V_aug)
  NULL
}

mean_or_na <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}


est_matrix <- function(res, method, T, value = "error", params = param_names()) {
  sub <- res[res$method == method & as.numeric(res$T) == as.numeric(T), , drop = FALSE]
  if (!nrow(sub)) return(matrix(numeric(0), 0, length(params), dimnames = list(NULL, params)))
  reps <- sort(unique(sub$rep))
  M <- matrix(NA_real_, length(reps), length(params), dimnames = list(reps, params))
  for (i in seq_along(reps)) {
    rr <- sub[sub$rep == reps[i], , drop = FALSE]
    for (p in params) {
      z <- rr[rr$param == p, value]
      if (length(z)) M[i, p] <- z[1]
    }
  }
  M[stats::complete.cases(M), , drop = FALSE]
}

summary_by_parameter <- function(res) {
  keys <- unique(res[, c("T", "method", "param")])
  rows <- list()
  for (i in seq_len(nrow(keys))) {
    key <- keys[i, ]
    sub <- res[as.numeric(res$T) == as.numeric(key$T) & res$method == key$method & res$param == key$param, , drop = FALSE]
    e <- sub$error
    rows[[length(rows) + 1L]] <- data.frame(
      scenario = bivar_spec()$name, T = as.numeric(key$T), method = key$method, param = key$param,
      group = unname(param_group(key$param)),
      n = sum(is.finite(e)), true = unique(sub$true)[1],
      bias = mean(e, na.rm = TRUE), sd = sd(e, na.rm = TRUE), rmse = sqrt(mean(e^2, na.rm = TRUE)),
      scaled_bias = sqrt(as.numeric(key$T)) * mean(e, na.rm = TRUE),
      scaled_sd = sqrt(as.numeric(key$T)) * sd(e, na.rm = TRUE),
      scaled_rmse = sqrt(as.numeric(key$T)) * sqrt(mean(e^2, na.rm = TRUE)),
      conv_rate = mean(sub$converged, na.rm = TRUE),
      boundary_rate = mean(sub$boundary_hit, na.rm = TRUE),
      at_lower_rate = mean(sub$at_lower, na.rm = TRUE), at_upper_rate = mean(sub$at_upper, na.rm = TRUE),
      info_pd_rate = mean(sub$info_pd, na.rm = TRUE),
      cover95 = mean(sub$cover95, na.rm = TRUE),
      ridge_active_rate = if ("ridge_active" %in% names(sub)) mean_or_na(as.logical(sub$ridge_active)) else NA_real_,
      mean_stage1_moment_norm = if ("stage1_moment_norm" %in% names(sub)) mean_or_na(sub$stage1_moment_norm) else NA_real_,
      mean_weighted_moment_norm = if ("weighted_moment_norm" %in% names(sub)) mean_or_na(sub$weighted_moment_norm) else NA_real_,
      mean_omega_cond = if ("omega_cond" %in% names(sub)) mean_or_na(sub$omega_cond) else NA_real_,
      mean_omega_min_eig = if ("omega_min_eig" %in% names(sub)) mean_or_na(sub$omega_min_eig) else NA_real_,
      mean_omega_n_eval = if ("omega_n_eval" %in% names(sub)) mean_or_na(sub$omega_n_eval) else NA_real_,
      mean_omega_n_intervals = if ("omega_n_intervals" %in% names(sub)) mean_or_na(sub$omega_n_intervals) else NA_real_,
      mean_aug_quad_hmax = if ("aug_quad_hmax" %in% names(sub)) mean_or_na(sub$aug_quad_hmax) else NA_real_,
      mean_aug_degree = if ("aug_degree" %in% names(sub)) mean_or_na(sub$aug_degree) else NA_real_,
      mean_q_aug = if ("q_aug" %in% names(sub)) mean_or_na(sub$q_aug) else NA_real_,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

standardized_summary <- function(res, pop) {
  ss <- summary_by_parameter(res)
  ss$target_sd_rootT <- NA_real_; ss$scaled_rmse_over_target <- NA_real_; ss$scaled_sd_over_target <- NA_real_
  for (i in seq_len(nrow(ss))) {
    V <- population_vcov(pop, ss$method[i])
    if (is.null(V)) next
    j <- match(ss$param[i], param_names())
    target <- sqrt(max(V[j, j], 0))
    ss$target_sd_rootT[i] <- target
    ss$scaled_rmse_over_target[i] <- ss$scaled_rmse[i] / target
    ss$scaled_sd_over_target[i] <- ss$scaled_sd[i] / target
  }
  ss
}


asymptotic_ci_long <- function(res, pop, level = 0.95) {
  level <- as.numeric(level)
  if (!is.finite(level) || level <= 0 || level >= 1) level <- 0.95
  zcrit <- qnorm(1 - (1 - level) / 2)
  params <- param_names()
  keys <- unique(res[, c("T", "method")])
  rows <- list()
  for (i in seq_len(nrow(keys))) {
    key <- keys[i, ]
    T0 <- as.numeric(key$T)
    method0 <- as.character(key$method)
    V <- population_vcov(pop, method0)
    if (is.null(V) || !all(dim(V) == c(length(params), length(params)))) next
    sd_rootT <- sqrt(pmax(diag(V), 0))
    names(sd_rootT) <- params
    sub <- res[as.numeric(res$T) == T0 & res$method == method0 & res$param %in% params, , drop = FALSE]
    if (!nrow(sub)) next
    j <- match(sub$param, params)
    target <- sd_rootT[j]
    asym_se <- target / sqrt(T0)
    half_width <- zcrit * asym_se
    width <- 2 * half_width
    cover <- abs(as.numeric(sub$error)) <= half_width
    cover[!is.finite(half_width)] <- NA
    standardized_error <- sqrt(T0) * as.numeric(sub$error) / target
    standardized_error[!is.finite(standardized_error)] <- NA_real_
    rows[[length(rows) + 1L]] <- data.frame(
      scenario = if ("scenario" %in% names(sub)) sub$scenario else bivar_spec()$name,
      T = T0,
      rep = sub$rep,
      method = method0,
      param = sub$param,
      group = unname(param_group(sub$param)),
      estimate = as.numeric(sub$estimate),
      true = as.numeric(sub$true),
      error = as.numeric(sub$error),
      ci_level = level,
      zcrit = zcrit,
      target_sd_rootT = as.numeric(target),
      asymptotic_se = as.numeric(asym_se),
      ci_half_width = as.numeric(half_width),
      ci_width = as.numeric(width),
      ci_lower = as.numeric(sub$estimate) - as.numeric(half_width),
      ci_upper = as.numeric(sub$estimate) + as.numeric(half_width),
      asymptotic_cover = as.logical(cover),
      standardized_error = as.numeric(standardized_error),
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

asymptotic_ci_by_parameter <- function(res, pop, level = 0.95) {
  ci <- asymptotic_ci_long(res, pop, level = level)
  if (!nrow(ci)) return(data.frame())
  keys <- unique(ci[, c("T", "method", "param")])
  rows <- list()
  for (i in seq_len(nrow(keys))) {
    key <- keys[i, ]
    sub <- ci[as.numeric(ci$T) == as.numeric(key$T) & ci$method == key$method & ci$param == key$param, , drop = FALSE]
    cov <- as.numeric(sub$asymptotic_cover)
    cov <- cov[is.finite(cov)]
    n <- length(cov)
    coverage <- if (n) mean(cov) else NA_real_
    mcse <- if (n && is.finite(coverage)) sqrt(coverage * (1 - coverage) / n) else NA_real_
    rows[[length(rows) + 1L]] <- data.frame(
      scenario = bivar_spec()$name,
      T = as.numeric(key$T),
      method = key$method,
      param = key$param,
      group = unname(param_group(key$param)),
      ci_level = unique(sub$ci_level)[1],
      zcrit = unique(sub$zcrit)[1],
      n = n,
      coverage = coverage,
      coverage_mcse = mcse,
      coverage_error = coverage - unique(sub$ci_level)[1],
      target_sd_rootT = unique(sub$target_sd_rootT)[1],
      asymptotic_se = unique(sub$asymptotic_se)[1],
      mean_ci_width = mean(sub$ci_width, na.rm = TRUE),
      median_ci_width = median(sub$ci_width, na.rm = TRUE),
      ci_half_width = unique(sub$ci_half_width)[1],
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$T, match(out$method, method_order()), match(out$param, param_names())), ]
}

asymptotic_ci_overall <- function(res, pop, level = 0.95) {
  byp <- asymptotic_ci_by_parameter(res, pop, level = level)
  if (!nrow(byp)) return(data.frame())
  byp$width_ratio_to_mle <- NA_real_
  for (T0 in sort(unique(byp$T))) {
    mle <- byp[byp$T == T0 & byp$method == "MLE", c("param", "mean_ci_width"), drop = FALSE]
    if (!nrow(mle)) next
    names(mle)[2] <- "mle_width"
    idx <- byp$T == T0
    tmp <- merge(byp[idx, c("method", "param", "mean_ci_width"), drop = FALSE], mle, by = "param", all.x = TRUE, sort = FALSE)
    tmp$ratio <- tmp$mean_ci_width / tmp$mle_width
    byp$width_ratio_to_mle[idx] <- tmp$ratio[match(paste(byp$method[idx], byp$param[idx]), paste(tmp$method, tmp$param))]
  }
  keys <- unique(byp[, c("T", "method")])
  rows <- list()
  for (i in seq_len(nrow(keys))) {
    key <- keys[i, ]
    sub <- byp[byp$T == key$T & byp$method == key$method, , drop = FALSE]
    rows[[length(rows) + 1L]] <- data.frame(
      scenario = bivar_spec()$name,
      T = as.numeric(key$T),
      method = key$method,
      ci_level = unique(sub$ci_level)[1],
      zcrit = unique(sub$zcrit)[1],
      n_parameters = nrow(sub),
      n_rep_parameter_cells = sum(sub$n, na.rm = TRUE),
      coverage = mean(sub$coverage, na.rm = TRUE),
      coverage_min = min(sub$coverage, na.rm = TRUE),
      coverage_max = max(sub$coverage, na.rm = TRUE),
      mean_coverage_mcse = mean(sub$coverage_mcse, na.rm = TRUE),
      mean_ci_width = mean(sub$mean_ci_width, na.rm = TRUE),
      median_ci_width = median(sub$median_ci_width, na.rm = TRUE),
      mean_width_ratio_to_mle = mean(sub$width_ratio_to_mle, na.rm = TRUE),
      median_width_ratio_to_mle = median(sub$width_ratio_to_mle, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$T, match(out$method, method_order())), ]
}

lan_quadratic_quantiles <- function(res, pop, probs = c(0.50, 0.90, 0.95, 0.99)) {
  keys <- unique(res[, c("T", "method")])
  rows <- list()
  for (i in seq_len(nrow(keys))) {
    key <- keys[i, ]
    V <- population_vcov(pop, key$method)
    if (is.null(V)) next
    Vinv <- safe_solve(V)
    E <- est_matrix(res, key$method, as.numeric(key$T), "error")
    if (nrow(E) < 5L) next
    Q <- as.numeric(as.numeric(key$T) * rowSums((E %*% Vinv) * E))
    qs <- quantile(Q, probs = probs, names = FALSE, na.rm = TRUE)
    rows[[length(rows) + 1L]] <- data.frame(
      scenario = bivar_spec()$name, T = as.numeric(key$T), method = key$method, prob = probs,
      q_empirical = qs, q_chisq = qchisq(probs, df = ncol(E)),
      ratio_to_chisq = qs / qchisq(probs, df = ncol(E)),
      n = sum(is.finite(Q)), stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

rate_slope_table <- function(res) {
  ss <- summary_by_parameter(res)
  keys <- unique(ss[, c("method", "param")])
  rows <- list()
  for (i in seq_len(nrow(keys))) {
    sub <- ss[ss$method == keys$method[i] & ss$param == keys$param[i] & ss$rmse > 0 & is.finite(ss$rmse), ]
    if (nrow(sub) < 3L) next
    fit <- lm(log(rmse) ~ log(T), data = sub)
    rows[[length(rows) + 1L]] <- data.frame(
      scenario = bivar_spec()$name, method = keys$method[i], param = keys$param[i],
      group = unname(param_group(keys$param[i])), n_T = nrow(sub),
      slope_log_rmse_on_log_T = unname(coef(fit)[2]),
      intercept = unname(coef(fit)[1]), r2 = summary(fit)$r.squared,
      target_slope = -0.5, stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

population_summary_table <- function(pop, B_boot = 1000, seed = 1, mat_long = NULL) {
  stats <- c("trace_ratio_ls_to_mle", "logdet_ratio_ls_to_mle",
             "eig_min_ls_to_mle", "eig_median_ls_to_mle", "eig_max_ls_to_mle")
  est <- c(pop$trace_ratio_ls, pop$logdet_ratio_ls,
           min(pop$eig_cov_ls_to_mle), median(pop$eig_cov_ls_to_mle), max(pop$eig_cov_ls_to_mle))
  if (!is.null(pop$V_aug)) {
    stats <- c(stats, "trace_ratio_aug_to_mle", "logdet_ratio_aug_to_mle",
               "eig_min_aug_to_mle", "eig_median_aug_to_mle", "eig_max_aug_to_mle")
    est <- c(est, pop$trace_ratio_aug, pop$logdet_ratio_aug,
             min(pop$eig_cov_aug_to_mle), median(pop$eig_cov_aug_to_mle), max(pop$eig_cov_aug_to_mle))
  }
  base <- data.frame(scenario = bivar_spec()$name, statistic = stats,
                     estimate = est, lo = NA_real_, hi = NA_real_, stringsAsFactors = FALSE)
  if (!is.null(mat_long) && length(unique(mat_long$replicate)) >= 3L && B_boot > 0L) {
    set.seed(seed)
    reps <- sort(unique(mat_long$replicate))
    boot <- replicate(B_boot, {
      rr <- sample(reps, replace = TRUE)
      xx <- do.call(rbind, lapply(seq_along(rr), function(j) {
        z <- mat_long[mat_long$replicate == rr[j], ]; z$replicate <- j; z
      }))
      pp <- build_population_list(xx)
      out <- c(trace_ratio_ls_to_mle = pp$trace_ratio_ls,
               logdet_ratio_ls_to_mle = pp$logdet_ratio_ls,
               eig_min_ls_to_mle = min(pp$eig_cov_ls_to_mle),
               eig_median_ls_to_mle = median(pp$eig_cov_ls_to_mle),
               eig_max_ls_to_mle = max(pp$eig_cov_ls_to_mle))
      if (!is.null(pp$V_aug)) {
        out <- c(out,
                 trace_ratio_aug_to_mle = pp$trace_ratio_aug,
                 logdet_ratio_aug_to_mle = pp$logdet_ratio_aug,
                 eig_min_aug_to_mle = min(pp$eig_cov_aug_to_mle),
                 eig_median_aug_to_mle = median(pp$eig_cov_aug_to_mle),
                 eig_max_aug_to_mle = max(pp$eig_cov_aug_to_mle))
      }
      out
    })
    ci <- t(apply(boot, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE))
    base$lo <- ci[base$statistic, 1]
    base$hi <- ci[base$statistic, 2]
  }
  base
}

population_parameter_table <- function(pop) {
  params <- param_names()
  V0 <- pop$I_inv
  V1 <- pop$V_ls
  V2 <- if (!is.null(pop$V_aug)) pop$V_aug else matrix(NA_real_, length(params), length(params))
  data.frame(
    scenario = bivar_spec()$name,
    param = params,
    group = unname(param_group(params)),
    sd_mle_rootT = sqrt(pmax(diag(V0), 0)),
    sd_gmm_rootT = sqrt(pmax(diag(V1), 0)),
    se_inflation_gmm_to_mle = sqrt(pmax(diag(V1), 0) / pmax(diag(V0), .Machine$double.eps)),
    var_inflation_gmm_to_mle = pmax(diag(V1), 0) / pmax(diag(V0), .Machine$double.eps),
    sd_aug_rootT = sqrt(pmax(diag(V2), 0)),
    se_inflation_aug_to_mle = sqrt(pmax(diag(V2), 0) / pmax(diag(V0), .Machine$double.eps)),
    var_inflation_aug_to_mle = pmax(diag(V2), 0) / pmax(diag(V0), .Machine$double.eps),
    stringsAsFactors = FALSE
  )
}

population_parameter_long_table <- function(pop) {
  wide <- population_parameter_table(pop)
  rbind(
    data.frame(scenario = wide$scenario, method = "GMM_Dtheta", param = wide$param, group = wide$group,
               sd_rootT = wide$sd_gmm_rootT,
               se_inflation_to_mle = wide$se_inflation_gmm_to_mle,
               var_inflation_to_mle = wide$var_inflation_gmm_to_mle,
               stringsAsFactors = FALSE),
    data.frame(scenario = wide$scenario, method = "GMM_aug", param = wide$param, group = wide$group,
               sd_rootT = wide$sd_aug_rootT,
               se_inflation_to_mle = wide$se_inflation_aug_to_mle,
               var_inflation_to_mle = wide$var_inflation_aug_to_mle,
               stringsAsFactors = FALSE)
  )
}

population_eigen_table <- function(pop) {
  out <- data.frame(
    scenario = bivar_spec()$name,
    method = "GMM_Dtheta",
    eig_index = seq_along(pop$eig_cov_ls_to_mle),
    covariance_eigenvalue = pop$eig_cov_ls_to_mle,
    information_eigenvalue = pop$eig_info_ls_to_mle,
    stringsAsFactors = FALSE
  )
  if (!is.null(pop$eig_cov_aug_to_mle)) {
    out <- rbind(out, data.frame(
      scenario = bivar_spec()$name,
      method = "GMM_aug",
      eig_index = seq_along(pop$eig_cov_aug_to_mle),
      covariance_eigenvalue = pop$eig_cov_aug_to_mle,
      information_eigenvalue = pop$eig_info_aug_to_mle,
      stringsAsFactors = FALSE
    ))
  }
  out
}

empirical_vs_population_covariance <- function(res, pop, B_boot = 500, seed = 1) {
  set.seed(seed)
  keys <- unique(res[, c("T", "method")])
  rows <- list()
  for (i in seq_len(nrow(keys))) {
    key <- keys[i, ]
    V <- population_vcov(pop, key$method)
    if (is.null(V)) next
    M <- est_matrix(res, key$method, as.numeric(key$T), "error")
    if (nrow(M) < 20L) next
    stat_fun <- function(ix) {
      Vhat <- symmetrize(cov(M[ix, , drop = FALSE]) * as.numeric(key$T))
      eig <- sort(eigen(symmetrize(sym_sqrt(V, inverse = TRUE) %*% Vhat %*% sym_sqrt(V, inverse = TRUE)), symmetric = TRUE, only.values = TRUE)$values)
      c(trace_ratio = sum(diag(Vhat)) / sum(diag(V)), eig_min = min(eig), eig_median = median(eig), eig_max = max(eig))
    }
    point <- tryCatch(stat_fun(seq_len(nrow(M))), error = function(e) rep(NA_real_, 4))
    boot <- replicate(B_boot, {
      ix <- sample.int(nrow(M), replace = TRUE)
      tryCatch(stat_fun(ix), error = function(e) rep(NA_real_, 4))
    })
    ci <- t(apply(boot, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE))
    rows[[length(rows) + 1L]] <- data.frame(
      scenario = bivar_spec()$name, T = as.numeric(key$T), method = key$method,
      statistic = names(point), estimate = as.numeric(point), lo = ci[, 1], hi = ci[, 2],
      n = nrow(M), stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

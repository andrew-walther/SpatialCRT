# ============================================================
# Script: methods_ext_participation.R
# Purpose: A3 partial participation + power and A4 unequal trained professionals
#   on the NC participant subsets (plan sections A3/A4), checkpointed per task.
# Author: Andrew Walther
# Created: 2026-10-20
# Dependencies: methods_ext_core.R, parallel, digest
# ============================================================
# Run from the project root (BLAS pinned to one thread; the runner refuses otherwise):
#   VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
#     Rscript application/code/methods_ext_participation.R <profile> <workers>
# Profiles: pilot (timing/smoke), a3 (production A3), a4 (production A4).
# Outputs: application/results/methods_ext_2026-10/participation/<profile>/
#   manifest.rds, checkpoints/<task>.rds (gitignored), run.log.
# Memory: one task holds at most 10,000 outcomes x ~70 fit cells of doubles (~6 MB),
# so 10 forked workers need well under 1 GB. On Longleaf, a SLURM array over the
# task list (one task per array index) would use the same checkpoints.
#
# DGP for participants P (queen, n = |P|), plan A3:
#   Y_P = (I - rho W_P)^{-1} [ tau Z + kappa tau e(Z) (x (1 - Z) control-only)
#                              + beta X_P + sigma o eps ],   eps ~ N(0, I)
# Fit: lean SAR lag on W_P with intercept, Z, unscaled e(Z) (x (1 - Z)), and X when
# beta = 1. Wald 95% interval and two-sided test of tau = 0.

.mp_dir <- local({
  for (i in rev(seq_len(sys.nframe()))) {
    p <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (!is.null(p)) return(normalizePath(dirname(p)))
  }
  f <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(dirname(sub("--file=", "", f[1]))) else stop("Run or source by path")
})
source(file.path(.mp_dir, "methods_ext_core.R"))

# Configuration ----
#' Return the approved A3/A4 settings, or the reduced timing pilot
#' @param profile "pilot", "a3" or "a4".
#' @return Configuration list (factors, replicate counts, effect grids).
#' @family methods_ext_participation
#' @examples
#' mp_config("pilot")
mp_config <- function(profile = c("pilot", "a3", "a4")) {
  profile <- match.arg(profile)
  base <- list(profile = profile, years = as.character(2018:2021), designs = 1:9,
               rho = c(0, 0.5), kappa = c(0.5, 0.8), regimes = c("both", "control_only"),
               S = 50L, J = 20L, R = 10L, J58 = 100L, R58 = 100L,
               det_R = 200L, det_R58 = 1000L, crit = stats::qnorm(0.975))
  if (profile == "a3") {
    base$N <- c(15, 20, 25, 29, 30, 35, 40, 45, 50, 58)
    base$mech <- list(random = base$N, incidence3x = c(20, 29))
    base$tau <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0)
    base$beta1 <- list(N = c(29, 58), tau = c(0, 0.4, 1.0))
    base$arms <- "unit"
  } else if (profile == "a4") {
    base$N <- c(20, 29, 40, 58)
    base$mech <- list(random = base$N)
    base$tau <- c(0, 0.1, 0.2, 0.3, 0.4)       # individual-SD effect d
    base$beta1 <- list(N = numeric(0), tau = numeric(0))
    base$arms <- c("equal_m", "unequal_m")
    base$m_bar <- 25; base$icc <- 0.05
  } else {
    base$years <- "2018"; base$N <- c(20, 58)
    base$mech <- list(random = c(20, 58))
    base$tau <- c(0, 1); base$kappa <- 0.8; base$rho <- 0.5
    base$beta1 <- list(N = numeric(0), tau = numeric(0)); base$arms <- "unit"
    base$S <- 10L; base$J <- 5L; base$R <- 10L; base$J58 <- 20L; base$R58 <- 20L
    base$det_R <- 50L; base$det_R58 <- 200L
  }
  base
}

#' Expand the task list: one task per (mechanism, year, N, design)
#' @param config From mp_config().
#' @return Data frame of tasks, heaviest (N = 58, then larger N) first.
#' @family methods_ext_participation
#' @examples
#' nrow(mp_tasks(mp_config("a3")))
mp_tasks <- function(config) {
  t <- do.call(rbind, lapply(names(config$mech), function(m)
    expand.grid(Mechanism = m, Year = config$years, N = config$mech[[m]],
                Design_ID = config$designs, stringsAsFactors = FALSE)))
  t$Task <- sprintf("%s_%s_N%02d_D%d", t$Mechanism, t$Year, t$N, t$Design_ID)
  t[order(-t$N, t$Design_ID %in% c(1, 2), t$Task), ]
}

#' Fit cells for one task: (beta, rho, regime, tau, kappa, arm)
#'
#' tau = 0 does not depend on kappa, so it is fitted once (kappa = NA) and
#' reported under both kappa values by the summaries (flag Shared_Tau0).
#' @param config Configuration.
#' @param mech Mechanism.
#' @param N Number of participants.
#' @return Data frame of cells.
#' @family methods_ext_participation
#' @examples
#' nrow(mp_cells(mp_config("a3"), "random", 29))
mp_cells <- function(config, mech, N) {
  one <- function(beta, taus) {
    pos <- expand.grid(tau = taus[taus > 0], kappa = config$kappa)
    tk <- rbind(if (0 %in% taus) data.frame(tau = 0, kappa = NA_real_), pos)
    do.call(rbind, lapply(config$arms, function(arm) do.call(rbind, lapply(config$rho, function(r)
      do.call(rbind, lapply(config$regimes, function(g)
        cbind(data.frame(Beta = beta, Rho = r, Regime = g, Arm = arm, stringsAsFactors = FALSE), tk)))))))
  }
  cells <- one(0, config$tau)
  if (mech == "random" && N %in% config$beta1$N) cells <- rbind(cells, one(1, config$beta1$tau))
  names(cells)[names(cells) == "tau"] <- "Tau"; names(cells)[names(cells) == "kappa"] <- "Kappa"
  rownames(cells) <- NULL
  cells
}

#' Per-college noise SD for an arm
#'
#' unit: sigma_c = 1 (A3; tau is the standardized cluster effect). equal_m:
#' sigma_c = sqrt(ICC + (1 - ICC)/m_bar). unequal_m: m_j = m_bar pop_j / mean(pop)
#' among participants, sigma_j = sqrt(ICC + (1 - ICC)/m_j) (A4).
#' @param arm Arm label.
#' @param pop Participants' populations (year's person-years).
#' @param config Configuration.
#' @return Numeric vector of SDs.
#' @family methods_ext_participation
#' @examples
#' mp_sigma("unit", c(1, 2), list())
mp_sigma <- function(arm, pop, config) {
  switch(arm, unit = rep(1, length(pop)),
         equal_m = rep(sqrt(config$icc + (1 - config$icc) / config$m_bar), length(pop)),
         unequal_m = { m <- config$m_bar * pop / mean(pop); sqrt(config$icc + (1 - config$icc) / m) },
         stop("Unknown arm"))
}

# One task ----
#' Simulate and fit every cell for one (mechanism, year, N, design) task
#'
#' For each subset s (shared across designs and, under random volunteering, years)
#' and allocation j, one standard-normal noise matrix (key excludes effects) feeds
#' every cell: Y = A_rho (tau z + kappa tau s_reg + beta X + sigma o eps). Per
#' (cell, subset, allocation) the counts and sums needed for power, size, coverage,
#' bias and MSE are kept; |t| at tau = 0 is kept per fit for size adjustment.
#' @param task One task row.
#' @param config Configuration.
#' @param setup Frozen setup.
#' @return List with cells, stats array [cell, subset, allocation, stat], t0 list,
#'   allocation diagnostics, subset diagnostics, fit count and elapsed seconds.
#' @family methods_ext_participation
#' @examples
#' # mp_task(mp_tasks(mp_config("pilot"))[1, ], mp_config("pilot"), me_load_setup())
mp_task <- function(task, config, setup) {
  t_start <- proc.time()[[3]]
  mech <- task$Mechanism; year <- task$Year; N <- task$N; id <- task$Design_ID
  cells <- mp_cells(config, mech, N)
  S <- if (N == 58) 1L else config$S
  a <- setup$annual[[year]]
  stat_names <- c("n_fit", "n_est", "n_ci", "n_rej", "n_cov", "sum_err", "sum_err2",
                  "n_alias", "n_bound", "n_error", "n_warn", "n_unident", "n_rej_adj")
  JR <- function(det) if (det) c(1L, if (N == 58) config$det_R58 else config$det_R) else
    if (N == 58) c(config$J58, config$R58) else c(config$J, config$R)
  stats <- NULL; tabs <- replicate(nrow(cells), list(), simplify = FALSE); alloc_rows <- list(); sub_rows <- list()
  n_fits <- 0
  for (s in seq_len(S)) {
    P <- me_subset(mech, N, s, year, setup)
    g <- me_geometry(P, setup)
    key <- if (mech == "random") list(mech, N, s) else list(mech, N, s, year)
    di <- me_design_inputs(g, year, setup, key)
    det <- me_deterministic(id, di); jr <- JR(det)
    if (is.null(stats)) {
      stats <- array(0, dim = c(nrow(cells), S, jr[1], length(stat_names)),
                     dimnames = list(NULL, NULL, NULL, stat_names))
    }
    if (dim(stats)[3] != jr[1]) stop("Allocation count differs across subsets in task ", task$Task)
    sub_rows[[s]] <- data.frame(Subset = s, N_Participants = g$n, Islands = g$islands,
      Components = g$components, Edges = g$n_edges, Engine_Error = g$engine_error,
      Listw_Warnings = paste(unique(g$listw_warnings), collapse = " | "), Guard_Bypass = isTRUE(g$guard_bypass),
      Min_Region_Size = min(tabulate(di$regions)), Regions_Below3 = sum(tabulate(di$regions) < 3),
      N_Regions = length(unique(di$regions)), N_Blocks = length(unique(di$blocks)),
      Min_Block_Size = min(table(di$blocks)), Deterministic = det, J = jr[1], R = jr[2],
      Participant_Deaths = sum(a$deaths[P]), Participant_PY = sum(a$person_years[P]),
      Mean_Participant_Rank = mean(a$rank[P]), stringsAsFactors = FALSE)
    pop <- a$person_years[P]
    sig <- lapply(config$arms, mp_sigma, pop = pop, config = config); names(sig) <- config$arms
    XP <- a$X[P]                                  # college's own rank-scaled burden (beta = 1)
    A <- lapply(config$rho, function(r) solve(diag(g$n) - r * g$WP)); names(A) <- config$rho
    for (j in seq_len(jr[1])) {
      do.call(me_seed, c(list("pp-Z"), key, list(year, id, j)))
      z <- me_draw_design(id, di)
      e <- me_exposure(z, g)
      alloc_rows[[length(alloc_rows) + 1L]] <- data.frame(Subset = s, Allocation = j,
        Treated = sum(z), Death_Share = sum(a$deaths[P] * z) / sum(a$deaths[P]),
        PY_Share = sum(a$person_years[P] * z) / sum(a$person_years[P]),
        Islands_Treated = sum(z[rowSums(g$AP) == 0]))
      do.call(me_seed, c(list("pp-eps"), key, list(year, id, j)))
      E <- matrix(stats::rnorm(g$n * jr[2]), nrow = g$n)
      for (reg in config$regimes) {
        sp <- if (reg == "control_only") e * (1 - z) else e
        xm0 <- cbind("(Intercept)" = 1, Z = z, Spill = sp)
        for (b in unique(cells$Beta)) {
          xm <- if (b == 1) cbind(xm0, X = XP) else xm0
          nuis <- xm[, colnames(xm) != "Z", drop = FALSE]
          ident <- sum(stats::lm.fit(nuis, z)$residuals^2) > 1e-10
          for (r in config$rho) {
            Ar <- A[[as.character(r)]]
            AZ <- as.vector(Ar %*% z); AS <- as.vector(Ar %*% sp); AX <- as.vector(Ar %*% (b * XP))
            for (arm in config$arms) {
              AE <- Ar %*% (sig[[arm]] * E)
              idx <- which(cells$Beta == b & cells$Rho == r & cells$Regime == reg & cells$Arm == arm)
              for (ci in idx) {
                tau <- cells$Tau[ci]; kap <- if (is.na(cells$Kappa[ci])) 0 else cells$Kappa[ci]
                mu <- tau * AZ + kap * tau * AS + AX
                res <- matrix(c(NA, NA, 0, 0, 1, 0), jr[2], 6, byrow = TRUE)   # no engine = error
                if (!is.null(g$engine)) for (k in seq_len(jr[2])) res[k, ] <- me_fit(mu + AE[, k], xm, g$engine)
                n_fits <- n_fits + jr[2]
                if (!ident) res[, 1:2] <- NA_real_
                err <- res[, 1] - tau; se <- res[, 2]
                ok <- is.finite(err); okci <- ok & is.finite(se) & se > 0
                tt <- abs(res[, 1] / se)
                tabs[[ci]][[length(tabs[[ci]]) + 1L]] <- list(s = s, j = j, t = ifelse(okci, tt, NA_real_))
                stats[ci, s, j, ] <- c(jr[2], sum(ok), sum(okci), sum(tt[okci] > config$crit),
                  sum(abs(err[okci]) <= config$crit * se[okci]), sum(err[ok]), sum(err[ok]^2),
                  sum(res[, 3] == 1), sum(res[, 4] == 1), sum(res[, 5] == 1),
                  sum(res[, 6] == 1), if (ident) 0 else jr[2], NA)
              }
            }
          }
        }
      }
    }
  }
  # Size-adjusted rejection: critical value = 95th percentile of |t| at tau = 0 for
  # the same (beta, rho, regime, arm); kept per fit so it can be recomputed later.
  t0 <- list()
  for (ci in which(cells$Tau == 0)) {
    grp <- which(cells$Beta == cells$Beta[ci] & cells$Rho == cells$Rho[ci] &
                 cells$Regime == cells$Regime[ci] & cells$Arm == cells$Arm[ci])
    tv <- unlist(lapply(tabs[[ci]], `[[`, "t"))
    crit <- if (any(is.finite(tv))) stats::quantile(tv, 0.95, na.rm = TRUE, type = 7) else NA_real_
    t0[[as.character(ci)]] <- list(cell = ci, crit = unname(crit), t = tv)
    for (cj in grp) for (rec in tabs[[cj]])
      stats[cj, rec$s, rec$j, "n_rej_adj"] <- sum(rec$t > crit, na.rm = TRUE)
  }
  list(task = task, cells = cells, stats = stats, t0 = t0,
       allocations = do.call(rbind, alloc_rows), subsets = do.call(rbind, sub_rows),
       n_fits = n_fits, elapsed = proc.time()[[3]] - t_start,
       config_hash = digest::digest(config))
}

# Driver ----
#' Run all tasks of a profile with forked workers and per-task checkpoints
#'
#' Finished tasks are skipped on restart (the manifest must match). Each worker
#' error is caught, logged and reported at the end; the run never stops silently.
#' @param profile Profile name.
#' @param workers Forked workers (10 locally).
#' @return Invisibly, a data frame of task status.
#' @family methods_ext_participation
#' @examples
#' # mp_run("pilot", 2)
mp_run <- function(profile = "pilot", workers = 2L) {
  me_check_threads()
  config <- mp_config(profile); setup <- me_load_setup()
  out <- file.path(.me_root, "participation", profile)
  me_manifest(out, config, c(file.path(.mp_dir, "methods_ext_core.R"), file.path(.mp_dir, "methods_ext_participation.R")))
  cp <- file.path(out, "checkpoints"); dir.create(cp, recursive = TRUE, showWarnings = FALSE)
  tasks <- mp_tasks(config)
  todo <- which(!file.exists(file.path(cp, paste0(tasks$Task, ".rds"))))
  cat(sprintf("[%s] %s: %d tasks, %d to run, %d workers\n", format(Sys.time()), profile, nrow(tasks), length(todo), workers))
  status <- parallel::mclapply(todo, function(i) {
    r <- tryCatch({
      res <- mp_task(tasks[i, ], config, setup)
      saveRDS(res, file.path(cp, paste0(tasks$Task[i], ".rds")))
      sprintf("OK %.0fs %d fits (%.3f ms/fit)", res$elapsed, res$n_fits, 1000 * res$elapsed / res$n_fits)
    }, error = function(e) paste("ERROR:", conditionMessage(e)))
    cat(sprintf("[%s] %s %s\n", format(Sys.time(), "%H:%M:%S"), tasks$Task[i], r))
    r
  }, mc.cores = workers, mc.preschedule = FALSE, mc.set.seed = FALSE)
  st <- data.frame(Task = tasks$Task[todo], Status = vapply(status, function(x)
    if (inherits(x, "try-error")) paste("WORKER FAILURE:", as.character(x)) else as.character(x), character(1)))
  write.csv(st, file.path(out, sprintf("status_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S"))), row.names = FALSE)
  bad <- !grepl("^OK", st$Status)
  cat(sprintf("[%s] done: %d OK, %d failed\n", format(Sys.time()), sum(!bad), sum(bad)))
  if (any(bad)) cat("FAILED TASKS:\n", paste(st$Task[bad], st$Status[bad], collapse = "\n"), "\n")
  invisible(st)
}

if (!exists("methods_ext_define_only") || !isTRUE(methods_ext_define_only)) {
  args <- commandArgs(TRUE)
  mp_run(if (length(args)) args[1] else "pilot", if (length(args) > 1) as.integer(args[2]) else 2L)
}

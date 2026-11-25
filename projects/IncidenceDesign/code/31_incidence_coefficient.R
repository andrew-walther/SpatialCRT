# ============================================================
# Script: 31_incidence_coefficient.R
# Purpose: Measure how the incidence coefficient affects estimation on the production grid.
# Author: Andrew Walther
# Created: 2026-11-25
# Dependencies: spdep, spatialreg, dplyr, digest, parallel (via 05)
# ============================================================
#
# From the IncidenceDesign root, use smoke, pilot or full:
#   VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
#     Rscript code/31_incidence_coefficient.R pilot 3
# The optional worker count defaults to 10. Each fork runs one production unit.
# Document notation: β multiplies Z (code tau), τ multiplies X (code beta).
# Output goes under results/incidence_coef_2026-11/<profile>/.

sim_define_only <- TRUE
sim_script_dir <- normalizePath("code", mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)
profile <- if (sys.nframe() == 0L && length(args)) args[[1L]] else "smoke"
stopifnot(profile %in% c("smoke", "pilot", "full"))
sim_profile <- if (profile == "smoke") "test" else profile
n_cores <- if (sys.nframe() == 0L && length(args) >= 2L) {
  stopifnot(length(args) == 2L, grepl("^([1-9]|10)$", args[[2L]]))
  as.integer(args[[2L]])
} else 10L
source("code/05_run_simulation.R")

incidence_values <- c(1, 0, 0.25, 0.5, 2)
scenario_keys <- c("Incidence_Mode", "Rho_Incidence", "Neighbor_Type", "Design",
                   "Rho", "Gamma", "Spillover_Type", "True_Tau", "Estimator")
production_tag <- "20260927_000502"

# Fit and summarize both coefficients ----

#' Fit one production model and keep the named Z and X coefficients
#' @param y Outcome vector, one entry per cluster.
#' @param xm Production model matrix with named Z and X columns.
#' @param setup Spatial lag setup from 04.
#' @param listw,engine Production fitter arguments. This study uses the lean engine.
#' @return List with coefficient estimates and SEs, residual df, and fit diagnostics.
#' @examples
#' # fit_incidence_model(y, cbind("(Intercept)" = 1, Z = Z, X = X), setup)
fit_incidence_model <- function(y, xm, setup, listw = NULL, engine = "lean") {
  stopifnot(engine == "lean")
  warns <- character()
  failed <- FALSE
  fit <- withCallingHandlers(
    tryCatch(fit_sar_lag(y, xm, setup), error = function(e) {
      failed <<- TRUE
      warns <<- c(warns, paste("ERROR:", conditionMessage(e)))
      NULL
    }),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  # Named subsetting returns NA when the fitter drops a coefficient.
  b <- stats::setNames(rep(NA_real_, 2L), c("Z", "X"))
  se <- b
  df <- NA_integer_
  if (!failed) {
    b[] <- fit$coefficients[names(b)]
    se[] <- fit$rest.se[names(se)]
    df <- as.integer(length(y) - length(fit$coefficients) - 1L)
  }
  list(tau = b[["Z"]], se = se[["Z"]], x_hat = b[["X"]], x_se = se[["X"]],
       df = df, warns = warns, failed = failed)
}

# Keep 04's model matrices and fit order. Redirect only this copy's per-model call.
fit_incidence_models <- local({
  fit_one_lag_model <- fit_incidence_model
  wrapper <- fit_tau_models
  environment(wrapper) <- environment()
  wrapper
})

#' Summarize one coefficient using production Wald and corrected t intervals
#' @param estimate,se Numeric vectors of estimates and fitted standard errors.
#' @param df Integer vector of per-fit degrees of freedom after dropping aliases.
#' @param surface,n_treated Surface index and treated count for each draw.
#' @param truth True coefficient.
#' @param two_sided If TRUE, count rejection when either bound excludes zero.
#' @return Scenario and surface data frames with production and corrected metrics.
#' @examples
#' # summarize_coefficient(c(0.8, 1.2), c(0.2, 0.2), c(95L, 95L), 1:2,
#' #                       c(50, 50), 1, TRUE)
summarize_coefficient <- function(estimate, se, df, surface, n_treated, truth, two_sided) {
  out <- summarize_fits(estimate, se, surface, n_treated, truth)
  valid <- is.finite(estimate) & is.finite(se)
  n <- length(surface)
  stopifnot(length(df) == n, all(df[valid] > 0))
  half <- qt(0.975, df) * sqrt(N / df) * se
  lo <- estimate - half
  hi <- estimate + half
  covered <- lo <= truth & hi >= truth
  rejected <- lo > 0 | (two_sided & hi < 0)
  out$scen$N_Attempted <- n
  out$scen$N_Missing_SE <- sum(is.finite(estimate) & !is.finite(se))
  out$scen$N_Missing_Estimate <- sum(!is.finite(estimate))
  out$scen$Mean_Model_SE <- if (any(valid)) mean(se[valid]) else NA_real_
  out$scen$Coverage_Tdf <- if (any(valid)) mean(covered[valid]) else NA_real_
  out$scen$Power_Tdf <- if (any(valid)) mean(rejected[valid]) else NA_real_
  for (k in seq_len(nrow(out$surf))) {
    on_surface <- surface == out$surf$Surface[k]
    v <- valid & on_surface
    out$surf$N_Attempted[k] <- sum(on_surface)
    out$surf$N_Missing_SE[k] <- sum(on_surface & is.finite(estimate) & !is.finite(se))
    out$surf$N_Missing_Estimate[k] <- sum(on_surface & !is.finite(estimate))
    out$surf$Mean_Estimate[k] <- if (any(v)) mean(estimate[v]) else NA_real_
    out$surf$SD[k] <- if (sum(v) > 1L) sd(estimate[v]) else NA_real_
    out$surf$Mean_Model_SE[k] <- if (any(v)) mean(se[v]) else NA_real_
    out$surf$Coverage_Tdf[k] <- if (any(v)) mean(covered[v]) else NA_real_
    out$surf$Power_Tdf[k] <- if (any(v)) mean(rejected[v]) else NA_real_
  }
  ok <- out$surf$N_Valid > 0
  out$scen$SE_Coverage_Tdf <- sd(out$surf$Coverage_Tdf[ok]) / sqrt(sum(ok))
  out$scen$SE_Power_Tdf <- sd(out$surf$Power_Tdf[ok]) / sqrt(sum(ok))
  out
}

# Simulation unit ----

#' Run one production unit at one incidence coefficient
#' @param cfg_index,nb_type,rho,X Production run_unit() arguments.
#' @param incidence_coef True τ, called beta in 05.
#' @param baseline List of surface-major Z estimates at τ = 1 for this unit.
#' @return Scenario, surface and warning tables and, at τ = 1, paired Z estimates.
#' @examples
#' # run_incidence_unit(1, "queen", 0, generate_surfaces(inc_configs[[1]]), 1)
run_incidence_unit <- function(cfg_index, nb_type, rho, X, incidence_coef, baseline = NULL) {
  cfg <- inc_configs[[cfg_index]]
  spatial <- get_active_spatial(grid_obj, nb_type)
  W <- spatial$W
  setup <- lag_setups[[nb_type]]
  inv_mat <- solve(I_mat - rho * W)
  KJ <- n_surfaces * n_design_draw
  surface <- rep(seq_len(n_surfaces), each = n_design_draw)
  scen_rows <- surf_rows <- warn_rows <- z_reference <- list()

  for (gamma_val in gamma_vals) {
    for (spill_type in spill_types) {
      set_seed_key("eps", cfg$mode, cfg$rho_x, nb_type, rho, gamma_val, spill_type)
      E <- matrix(rnorm(N * KJ, 0, sigma), nrow = N, ncol = KJ)
      for (d_id in design_ids) {
        Zm <- draw_assignments(cfg, nb_type, rho, gamma_val, spill_type, d_id, X)
        WZ <- W %*% Zm
        Sm <- if (spill_type == "control_only") gamma_val * WZ * (1 - Zm) else gamma_val * WZ
        Xm <- X[, surface]
        n_treated <- colSums(Zm)
        rank_def <- any(vapply(seq_len(KJ), function(f) {
          qr(cbind(1, Zm[, f], WZ[, f]))$rank < 3
        }, logical(1)))
        # Correlation describes assignment alignment, including draws with unusable fits.
        x_sd <- vapply(seq_len(n_surfaces), function(k) sd(X[, k]), numeric(1))
        abs_cor <- vapply(seq_len(KJ), function(f) {
          if (n_treated[f] %in% c(0, N) || x_sd[surface[f]] == 0) return(NA_real_)
          abs(cor(Zm[, f], Xm[, f]))
        }, numeric(1))
        cor_k <- vapply(split(abs_cor, surface), mean, numeric(1), na.rm = TRUE)
        n_cor_k <- as.integer(rowsum(as.integer(is.finite(abs_cor)), surface))

        # Preserve 05's operation order exactly for Y = β(I - ρW)^(-1)Z + base.
        base <- inv_mat %*% (Sm + incidence_coef * Xm + E)
        IZ <- inv_mat %*% Zm
        for (true_tau in true_tau_vals) {
          Y <- true_tau * IZ + base
          fits <- lapply(seq_len(KJ), function(f) {
            fit_incidence_models(Y[, f], Zm[, f], Sm[, f], Xm[, f], setup,
                                 listw = spatial$listw, engine = engine)
          })
          for (est in c("oracle", "nonoracle")) {
            model <- lapply(fits, `[[`, est)
            z_hat <- vapply(model, `[[`, numeric(1), "tau")
            z_se <- vapply(model, `[[`, numeric(1), "se")
            x_hat <- vapply(model, `[[`, numeric(1), "x_hat")
            x_se <- vapply(model, `[[`, numeric(1), "x_se")
            df <- vapply(model, `[[`, integer(1), "df")
            warns <- lapply(model, `[[`, "warns")
            failed <- vapply(model, `[[`, logical(1), "failed")
            is_alias <- vapply(warns, function(w) any(startsWith(w, "Aliased variables found")),
                               logical(1))
            other <- lapply(warns, function(w) w[!startsWith(w, "Aliased variables found")])
            keys <- data.frame(
              Incidence_Mode = cfg$mode, Rho_Incidence = cfg$rho_x,
              Neighbor_Type = nb_type, Design = paste("Design", d_id),
              Rho = rho, Gamma = gamma_val, Spillover_Type = spill_type,
              True_Tau = true_tau, Estimator = est, Incidence_Coef = incidence_coef)
            z <- summarize_coefficient(z_hat, z_se, df, surface, n_treated, true_tau, FALSE)
            x <- summarize_coefficient(x_hat, x_se, df, surface, n_treated, incidence_coef, TRUE)
            # 05 also returns one-sided Wald rejection for X. Store only two-sided
            # rejection for X so it has one definition across result columns.
            x$scen[c("Power", "SE_Power")] <- NULL
            x$surf$Power <- NULL

            pair_key <- paste(gamma_val, spill_type, d_id, true_tau, est)
            z_pair <- replace(z_hat, !is.finite(z_se), NA_real_)
            if (incidence_coef == 1) {
              z_reference[[pair_key]] <- z_pair
              ref <- z_pair
            } else {
              ref <- baseline[[pair_key]]
              if (length(ref) != KJ) stop("Missing paired baseline for ", pair_key)
            }
            delta <- z_pair - ref
            delta[!is.finite(delta)] <- NA_real_
            delta_k <- split(delta, surface)
            z$scen$N_Paired <- sum(is.finite(delta))
            z$scen$Paired_MSD <- if (z$scen$N_Paired) mean(delta^2, na.rm = TRUE) else NA_real_
            z$surf$N_Paired <- vapply(delta_k, function(d) sum(is.finite(d)), integer(1))
            z$surf$Paired_MSD <- vapply(delta_k, function(d) {
              if (any(is.finite(d))) mean(d^2, na.rm = TRUE) else NA_real_
            }, numeric(1))
            flags <- data.frame(N_Aliased = sum(is_alias), N_Warn = sum(lengths(other)),
                                Z_WZ_rank_deficient = rank_def, N_Failed = sum(failed),
                                N_Warned = sum(lengths(warns) > 0L & !failed),
                                N_Other_Warned = sum(vapply(other, function(w) {
                                  any(!startsWith(w, "ERROR:"))
                                }, logical(1))),
                                Mean_Abs_Cor_ZX = if (any(is.finite(abs_cor))) {
                                  mean(abs_cor, na.rm = TRUE)
                                } else NA_real_, N_Cor_ZX = sum(n_cor_k))
            alignment <- data.frame(Mean_Abs_Cor_ZX = ifelse(n_cor_k > 0L, cor_k, NA_real_),
                                    N_Cor_ZX = n_cor_k, row.names = NULL)
            for (coef in c("Z", "X")) {
              sm <- if (coef == "Z") z else x
              scen_rows[[length(scen_rows) + 1L]] <- cbind(keys, Coefficient = coef, sm$scen, flags)
              surf_rows[[length(surf_rows) + 1L]] <- cbind(keys, Coefficient = coef, sm$surf, alignment)
            }
            tab <- table(unlist(warns))
            if (length(tab)) {
              warn_rows[[length(warn_rows) + 1L]] <- cbind(
                keys, Message = names(tab), Count = as.integer(tab))
            }
          }
        }
      }
    }
  }
  list(scen = bind_rows(scen_rows), surf = bind_rows(surf_rows),
       warn = bind_rows(warn_rows), baseline = z_reference)
}

# Keys and production reproduction ----

#' Match a complete set of rows by key and return actual indices in expected order
#' @param expected,actual Data frames containing the key columns.
#' @param keys Character vector of key column names.
#' @return Integer vector of actual row indices. Stops on missing, extra or duplicate keys.
#' @examples
#' match_grid(data.frame(id = 1:2), data.frame(id = 2:1), "id")
match_grid <- function(expected, actual, keys) {
  a <- expected[keys]
  b <- actual[keys]
  if (anyNA(a) || anyNA(b) || anyDuplicated(a) || anyDuplicated(b) || nrow(a) != nrow(b)) {
    stop("Invalid grid keys or row count. Expected ", nrow(a), ", got ", nrow(b))
  }
  idx <- dplyr::left_join(a, transform(b, Matched_Row = seq_len(nrow(b))), by = keys)$Matched_Row
  if (anyNA(idx)) stop("Missing grid keys: ", sum(is.na(idx)))
  idx
}

#' Compare the τ = 1 Z metrics to the dated production outputs
#' @param scen,surf Study tables at τ = 1 with both coefficients and estimators.
#' @param out_dir Profile's results directory for the reproduction report.
#' @return Invisibly TRUE on exact agreement. Writes a report before stopping on mismatch.
#' @examples
#' # check_reproduction(scenario_rows, surface_rows, "results/incidence_coef_2026-11/pilot")
check_reproduction <- function(scen, surf, out_dir) {
  tags <- c("MLE_tau_sweep", "MLEnonoracle_tau_sweep")
  lines <- sprintf("Production tag %s | profile %s | incidence coefficient = 1", production_tag, profile)
  ok <- TRUE
  for (table_name in c("scenarios", "surfaces")) {
    actual <- if (table_name == "scenarios") scen else surf
    actual <- actual[actual$Coefficient == "Z", ]
    keys <- if (table_name == "scenarios") scenario_keys else c(scenario_keys, "Surface")
    patterns <- c(scenarios = "sim_results_%s_combined_%s.rds", surfaces = "surface_results_%s_%s.rds")
    expected <- bind_rows(lapply(tags, function(tag) {
      readRDS(file.path("results/sim_data", sprintf(patterns[[table_name]], tag, production_tag)))
    }))
    if (profile == "pilot") {
      expected <- expected[expected$Rho %in% rho_vals & expected$Gamma %in% gamma_vals &
                             expected$True_Tau %in% true_tau_vals, ]
    }
    lines <- c(lines, sprintf("%s: %d production rows, %d rerun rows", table_name,
                              nrow(expected), nrow(actual)))
    idx <- tryCatch(match_grid(expected, actual, keys), error = identity)
    if (inherits(idx, "error")) {
      lines <- c(lines, conditionMessage(idx))
      ok <- FALSE
      next
    }
    metrics <- setdiff(names(expected), keys)
    equal <- vapply(metrics, function(m) {
      identical(expected[[m]], actual[[m]][idx], num.eq = FALSE, single.NA = FALSE)
    }, logical(1))
    lines <- c(lines, sprintf("  %-25s %s", metrics, ifelse(equal, "exact", "MISMATCH")))
    ok <- ok && all(equal)
  }
  lines <- c(lines, if (ok) "PASS: all production metrics match exactly." else "FAIL: production mismatch.")
  writeLines(lines, file.path(out_dir, "reproduction_check.txt"))
  cat(lines, sep = "\n")
  if (!ok) stop("τ = 1 reproduction failed. See reproduction_check.txt.")
  invisible(TRUE)
}

# Run and checkpoint ----

#' Run the selected profile and write verified scenario, surface and diagnostic tables
#' @return Invisibly, the profile output directory.
#' @examples
#' # Run from the command line with BLAS pinned as shown in the script header.
#' # run_incidence_study()
run_incidence_study <- function() {
  blas_vars <- c("VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS")
  if (!all(Sys.getenv(blas_vars) == "1")) stop("Set all three BLAS thread variables to 1 before R starts.")
  out_dir <- file.path("results/incidence_coef_2026-11", profile)
  cp_dir <- file.path(out_dir, "checkpoints")
  dir.create(cp_dir, recursive = TRUE, showWarnings = FALSE)
  started <- Sys.time()

  params <- list(profile = profile, engine = engine, grid_dim = grid_dim,
                 inc_configs = inc_configs, base_rate = base_rate, pop_per_cluster = pop_per_cluster,
                 pop_mode = pop_mode, rho_vals = rho_vals, gamma_vals = gamma_vals,
                 spill_types = spill_types, neighbor_types = neighbor_types,
                 true_tau_vals = true_tau_vals, incidence_values = incidence_values, sigma = sigma,
                 n_surfaces = n_surfaces, n_design_draw = n_design_draw, design_ids = design_ids)
  files <- c(code_files, "code/31_incidence_coefficient.R")
  manifest <- list(params = params,
                   code = setNames(vapply(files, digest::digest, character(1), file = TRUE),
                                   basename(files)),
                   packages = vapply(c("spatialreg", "spdep", "dplyr", "digest"), function(p) {
                     as.character(packageVersion(p))
                   }, character(1)),
                   R = R.version.string, blas = sessionInfo()$BLAS, lapack = extSoftVersion()["LAPACK"])
  manifest_file <- file.path(cp_dir, "manifest.rds")
  if (file.exists(manifest_file)) {
    if (!identical(readRDS(manifest_file), manifest)) {
      stop("Checkpoint manifest mismatch in ", cp_dir, ". Move the old checkpoint folder before rerunning.")
    }
  } else {
    if (length(list.files(cp_dir, pattern = "^unit_"))) stop("Checkpoints have no manifest in ", cp_dir)
    saveRDS(manifest, manifest_file)
  }

  units <- expand.grid(cfg_index = seq_along(inc_configs), nb_type = neighbor_types,
                       rho = rho_vals, stringsAsFactors = FALSE)
  surfaces <- lapply(inc_configs, generate_surfaces)
  expected <- bind_rows(lapply(inc_configs, function(cfg) {
    expand.grid(Incidence_Mode = cfg$mode, Rho_Incidence = cfg$rho_x,
                Neighbor_Type = neighbor_types, Design = paste("Design", design_ids),
                Rho = rho_vals, Gamma = gamma_vals, Spillover_Type = spill_types,
                True_Tau = true_tau_vals, Estimator = c("oracle", "nonoracle"),
                Coefficient = c("Z", "X"), stringsAsFactors = FALSE)
  }))
  expected_surfaces <- merge(expected, data.frame(Surface = seq_len(n_surfaces)))
  timings <- numeric(length(incidence_values))
  results <- vector("list", length(incidence_values))
  cat(sprintf("Profile %s | %d forked workers | %d units per incidence coefficient\n",
              profile, n_cores, nrow(units)))

  for (j in seq_along(incidence_values)) {
    value <- incidence_values[[j]]
    t0 <- Sys.time()
    cat(sprintf("Starting incidence coefficient %g\n", value))
    # Store only τ = 1 per-fit Z estimates, in the ignored baseline checkpoints.
    unit_results <- parallel::mclapply(seq_len(nrow(units)), function(i) {
      u <- units[i, ]
      suffix <- sprintf("cfg%02d_%s_rho%.2f.rds", u$cfg_index, u$nb_type, u$rho)
      cp <- file.path(cp_dir, sprintf("unit_tau%.2f_%s", value, suffix))
      if (file.exists(cp)) {
        out <- readRDS(cp)
      } else {
        baseline <- if (value != 1) {
          readRDS(file.path(cp_dir, paste0("unit_tau1.00_", suffix)))$baseline
        }
        out <- run_incidence_unit(u$cfg_index, u$nb_type, u$rho,
                                  surfaces[[u$cfg_index]], value, baseline)
        tmp <- tempfile(pattern = "pending_", tmpdir = cp_dir)
        saveRDS(out, tmp)
        if (!file.rename(tmp, cp)) stop("Could not publish checkpoint ", cp)
      }
      out$baseline <- NULL
      out
    }, mc.cores = n_cores, mc.preschedule = FALSE)
    bad <- vapply(unit_results, function(r) is.null(r) || inherits(r, "try-error"), logical(1))
    if (any(bad)) {
      msgs <- vapply(unit_results[bad], function(r) {
        if (is.null(r)) "worker returned no result" else as.character(r)
      }, character(1))
      stop("Failed production units: ", paste(which(bad), msgs, sep = " = ", collapse = "; "))
    }

    scen <- bind_rows(lapply(unit_results, `[[`, "scen"))
    surf <- bind_rows(lapply(unit_results, `[[`, "surf"))
    warn <- bind_rows(lapply(unit_results, `[[`, "warn"))
    match_grid(expected, scen, c(scenario_keys, "Coefficient"))
    match_grid(expected_surfaces, surf, c(scenario_keys, "Coefficient", "Surface"))
    stopifnot(all(scen$Incidence_Coef == value), all(surf$Incidence_Coef == value),
              all(scen$N_Attempted == n_surfaces * n_design_draw),
              all(surf$N_Attempted == n_design_draw),
              all(scen$N_Valid_Est + scen$N_Missing_Estimate + scen$N_Missing_SE == scen$N_Attempted),
              all(surf$N_Valid + surf$N_Missing_Estimate + surf$N_Missing_SE == surf$N_Attempted))
    if (value == 1 && profile != "smoke") check_reproduction(scen, surf, out_dir)
    z <- scen[scen$Coefficient == "Z", ]
    x <- scen[scen$Coefficient == "X", ]
    timings[[j]] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf(paste0("τ = %g | %d scenario-estimator rows | %d errors | %d aliased fits | ",
                       "%d fits with other warnings | unavailable Z/X: %d/%d | %.1f seconds\n"),
                value, nrow(z), sum(z$N_Failed), sum(z$N_Aliased), sum(z$N_Other_Warned),
                sum(z$N_Missing_Estimate + z$N_Missing_SE),
                sum(x$N_Missing_Estimate + x$N_Missing_SE), timings[[j]]))
    results[[j]] <- list(scen = scen, surf = surf, warn = warn)
  }

  scen <- bind_rows(lapply(results, `[[`, "scen"))
  surf <- bind_rows(lapply(results, `[[`, "surf"))
  warn <- bind_rows(lapply(results, `[[`, "warn"))
  if (!nrow(warn)) {
    warn <- cbind(scen[FALSE, c(scenario_keys, "Incidence_Coef")], Message = character(), Count = integer())
  }
  saveRDS(scen, file.path(out_dir, "scenarios.rds"))
  saveRDS(surf, file.path(out_dir, "surfaces.rds"))
  write.csv(warn, file.path(out_dir, "fit_messages.csv"), row.names = FALSE)
  saveRDS(list(manifest = manifest, n_cores = n_cores, session = sessionInfo(),
               timing = data.frame(Incidence_Coef = incidence_values, Wall_Seconds = timings)),
          file.path(out_dir, "run_info.rds"))
  cat(sprintf("Saved %d scenario rows and %d surface rows in %.1f min to %s\n",
              nrow(scen), nrow(surf), as.numeric(difftime(Sys.time(), started, units = "mins")), out_dir))
  invisible(out_dir)
}

if (sys.nframe() == 0L) run_incidence_study()

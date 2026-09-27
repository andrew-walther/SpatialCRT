# ============================================================
# Script: test_allocation_risk_summary.R
# Purpose: Verify fixed-unit aggregation, MC errors, SRS ratios and tail diagnostics.
# Author: Andrew Walther
# Created: 2026-09-27
# Dependencies: base R
# ============================================================
f <- grep("--file=", commandArgs(FALSE), value = TRUE)[1]
source(file.path(dirname(dirname(normalizePath(sub("--file=", "", f)))), "17_allocation_risk_summary.R"))
# Two fixed units; independent design streams. Deliberately unequal SRS means
# make the ratio-of-averages differ from the average of unit-specific ratios.
s <- data.frame(Incidence_Mode = "iid", Rho_Incidence = 0, Surface = rep(1:2, 2),
  Neighbor_Type = "queen", Rho = 0.5, Gamma = 0.8, Spillover_Type = "both",
  True_Tau = 1, Estimator = "oracle", Design = rep(c("Design 6", "Design 9"), each = 2),
  Design_Name = rep(c("Design 6: Balanced Quartiles", "Design 9: Simple Random Sampling"), each = 2),
  Complete = TRUE, Mean_MSE = c(1, 3, 2, 10), SE_Mean_MSE_Joint = c(3, 4, 6, 8),
  Variance_Corrected = c(-1, 3, 4, 4), Q90_Estimated = c(2, 4, 3, 7),
  Worst10_Mean_Estimated = c(3, 5, 4, 8), Sampled_Max_Estimated = c(4, 6, 5, 9),
  Coverage = c(.9, 1, .8, 1))
g <- "Spillover_Type"
x <- allocation_risk_compare(s, g); d <- x[x$Design == "Design 6", ]
stopifnot(d$Mean_Allocation_MSE == 2, d$MCSE_Mean_MSE == 2.5,
  d$Average_Corrected_Within_Variance == 1, d$RMS_Corrected_Within_SD == 1,
  d$N_Negative_Variance == 1, d$Average_Conditional_Q90 == 3,
  d$Average_Conditional_Worst10 == 4, d$Observed_Max_Sample_Dependent == 6,
  d$Mean_Allocation_MSE_Ratio_To_SRS == 1/3, d$Mean_Difference_To_SRS == -4,
  abs(d$MCSE_Difference_To_SRS - sqrt(2.5^2 + 5^2)) < 1e-12,
  d$Share_Fixed_Units_Below_SRS_Mean_MSE == 1,
  x$MCSE_Difference_To_SRS[x$Design == "Design 9"] == 0)
# Group isolation: huge losses in another regime cannot contaminate this result.
s2 <- s; s2$Spillover_Type <- "control_only"; s2$Mean_MSE <- s2$Mean_MSE * 100
z <- allocation_risk_compare(rbind(s, s2), g)
stopifnot(z$Mean_Allocation_MSE[z$Design == "Design 6" & z$Spillover_Type == "both"] == 2,
  z$Mean_Allocation_MSE[z$Design == "Design 6" & z$Spillover_Type == "control_only"] == 200)
# Incomplete and unmatched designs fail loudly rather than dropping rows.
bad <- s; bad$Complete[1] <- FALSE
stopifnot(inherits(try(allocation_risk_compare(bad, g), silent = TRUE), "try-error"),
  inherits(try(allocation_risk_compare(s[-1, ], g), silent = TRUE), "try-error"))
# Half-selected top draw reverses completely: selected=10, held-out=0.
a <- s[rep(1, 2), c(allocation_unit_keys, "Design")]
a$Frequency <- c(1, 1); a$MSE_Half1 <- c(10, 0); a$MSE_Half2 <- c(0, 10)
p <- s[1, ]; p$Q90_Half1 <- 9; p$Q90_Half2 <- 9
p$Worst10_Half1 <- 10; p$Worst10_Half2 <- 10; p$Tail_Half_Jaccard <- 0
p$Tail_Half_Rank_Cor <- -1; p$Tail_Boundary_Uncertain <- TRUE
v <- allocation_risk_precision(a, p)
stopifnot(v$Half_Selected_Worst10 == 10, v$Cross_Selected_Worst10 == 0,
          v$Selection_Noise_Optimism == 10)
# Unequal draw frequencies determine selection; never count unique caches equally.
a$Frequency <- c(9, 1); a$MSE_Half1 <- c(1, 10); a$MSE_Half2 <- c(2, 20)
v <- allocation_risk_precision(a, p)
stopifnot(v$Half_Selected_Worst10 == 15, v$Cross_Selected_Worst10 == 15)
# Replication comparison matches identities even when original rows are reordered.
s$SD_Corrected <- c(NA, sqrt(3), 2, 2); s$Variance_MC_Noise <- .5
s$Variance_Noise_Fraction <- .1; s$Tail_Half_Jaccard <- .5; s$Tail_Half_Rank_Cor <- .8
extended <- s[2:4, ]; extended$Mean_MSE <- extended$Mean_MSE * 2
v <- allocation_risk_extension(extended, s[4:1, ])
stopifnot(all(v$Mean_MSE_Ratio == 2), identical(v$Mean_MSE_Pilot, s$Mean_MSE[2:4]))
wrong <- extended; wrong$Surface[1] <- 99
stopifnot(inherits(try(allocation_risk_extension(wrong, s), silent = TRUE), "try-error"))
cat("Allocation-risk summary behavior checks passed\n")

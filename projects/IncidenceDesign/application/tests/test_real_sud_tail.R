# ============================================================
# Script: test_real_sud_tail.R
# Purpose: Check cross-selected risk diagnostics without post-selection reuse.
# Author: Andrew Walther
# Created: 2026-10-02
# Dependencies: existing real_sud tail helpers
# ============================================================
test_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
real_sud_tail_define_only <- TRUE
source(file.path(test_dir, "..", "code", "refine_real_sud_tail.R"))
# The halves disagree about the highest-risk allocation: selecting and evaluating
# in the same half gives 9; evaluation in its independent other half gives 1.
b <- list(allocations = data.frame(Frequency = c(5, 5), Complete = TRUE,
  MSE_Half1 = c(9, 1), MSE_Half2 = c(1, 9)))
x <- rs_cross_tail(b)
stopifnot(x$Tail_Draws == 1, x$Worst10_Select1_Evaluate2 == 1,
  x$Worst10_Select2_Evaluate1 == 1)
b$allocations$Complete[1] <- FALSE
stopifnot(is.na(rs_cross_tail(b)$Worst10_Select1_Evaluate2))
b$allocations <- data.frame(Frequency = 100, Complete = TRUE, MSE_Half1 = 2, MSE_Half2 = 3)
x <- rs_cross_tail(b)
stopifnot(x$Tail_Draws == 10, x$Worst10_Select1_Evaluate2 == 3,
  x$Worst10_Select2_Evaluate1 == 2)
cat("PASS: cross-selected tail values use independent halves, respect frequencies and withhold incomplete risks.\n")

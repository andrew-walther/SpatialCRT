# =============================================================================
# install_packages.R — First-time R package setup for Longleaf
#
# Run this ONCE in an interactive session before submitting jobs.
#
# Usage:
#   srun -p interact -n 1 --mem=8g -t 1:00:00 --pty /bin/bash
#   module add r/4.4.0
#   Rscript install_packages.R
# =============================================================================

cat("Installing R packages for SpatialCRT simulation...\n")
cat("R version:", R.version.string, "\n")
cat("Library paths:", .libPaths(), "\n\n")

# Required packages
packages <- c("spdep", "spatialreg", "dplyr", "digest")

# Install any that are missing
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, repos = "https://cloud.r-project.org")
  } else {
    cat(pkg, "already installed.\n")
  }
}

# Verify all packages load correctly
cat("\nVerifying package loading:\n")
all_ok <- TRUE
for (pkg in packages) {
  result <- tryCatch({
    library(pkg, character.only = TRUE)
    cat("  ", pkg, "- OK (version", as.character(packageVersion(pkg)), ")\n")
    TRUE
  }, error = function(e) {
    cat("  ", pkg, "- FAILED:", conditionMessage(e), "\n")
    FALSE
  })
  if (!result) all_ok <- FALSE
}

if (all_ok) {
  cat("\nAll packages installed and verified successfully.\n")
} else {
  cat("\nSome packages failed to load.\n")
  cat("If sf/spdep failed, try loading system libraries first:\n")
  cat("  module add gdal geos proj\n")
  cat("Then re-run this script.\n")
}

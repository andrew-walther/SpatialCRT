# ==============================================================================
# DIAGNOSTIC SANDBOX: Single Iteration Test
# Goal: Isolate the "spdep" crash without parallel overhead.
# ==============================================================================

library(tidyverse)
library(sf)
library(spdep)
library(spatialreg)

# 1. SETUP: Small 2x3 Grid (Simplest case that still allows blocking)
rows <- 2
cols <- 3
pts_per_cell <- 20 # Random variation allowed later
total_districts <- rows * cols
pts_total <- total_districts * pts_per_cell

cat("--- STEP 1: Grid Setup ---\n")
# Create Grid
grid_geom <- st_make_grid(
  st_as_sfc(st_bbox(c(xmin = 0, xmax = cols * 100, ymin = 0, ymax = rows * 100))),
  n = c(cols, rows), what = "polygons", square = TRUE
)
grid_sf <- st_sf(district = 1:total_districts, geometry = grid_geom)
nb_dist <- spdep::cell2nb(nrow = rows, ncol = cols, type = "rook")
mat_dist <- spdep::nb2mat(nb_dist, style = "B", zero.policy = TRUE)

cat("Grid created. District Adjacency Matrix dimensions:", dim(mat_dist), "\n")

# 2. DATA GENERATION: Rejection Sampling
cat("\n--- STEP 2: Point Generation (Rejection Sampling) ---\n")

valid_sample <- FALSE
points_df <- NULL

# Attempt to generate sample until no districts are empty
attempt <- 0
while(!valid_sample) {
  attempt <- attempt + 1
  d_samples <- sample(1:total_districts, pts_total, replace = TRUE)
  
  # Check coverage
  if(length(unique(d_samples)) == total_districts) {
    points_df <- data.frame(district = d_samples)
    points_df$id <- 1:nrow(points_df)
    valid_sample <- TRUE
    cat("Success! Generated valid distribution on attempt:", attempt, "\n")
  }
}

# Add dummy treatment (just for the test)
points_df$intervention <- sample(0:1, pts_total, replace = TRUE)
points_df$spillover <- 0 # Dummy value

# 3. WEIGHT MATRIX CONSTRUCTION (The Danger Zone)
cat("\n--- STEP 3: Constructing Point-Level Weights ---\n")

# Map district matrix to point matrix
d_idx <- points_df$district
W_points <- mat_dist[d_idx, d_idx]

# SANITIZATION CHECK 1
cat("Matrix created. Checking for NAs/NaNs in raw matrix...\n")
if(any(is.na(W_points))) cat("WARNING: NAs found in W_points!\n")
if(any(is.nan(W_points))) cat("WARNING: NaNs found in W_points!\n")

# Clean attributes
W_points <- as.matrix(W_points)
attr(W_points, "dimnames") <- NULL
storage.mode(W_points) <- "double"

cat("Matrix dimensions:", dim(W_points), "\n")

# 4. CONVERT TO LISTW (The Crash Site)
cat("\n--- STEP 4: Converting Matrix to listw ---\n")

# Test A: Binary Weights
cat("Attempting Binary (style='B') conversion...\n")
lw_bin <- tryCatch({
  mat2listw(W_points, style = "B", zero.policy = TRUE)
}, error = function(e) {
  cat("ERROR in Binary Conversion:", e$message, "\n")
  return(NULL)
})

if(!is.null(lw_bin)) {
  cat("Binary conversion successful.\n")
  
  # Check for islands (points with 0 neighbors)
  cards <- card(lw_bin$neighbours)
  num_islands <- sum(cards == 0)
  cat("Number of 'Islands' (points with 0 neighbors):", num_islands, "\n")
  
  # Test B: Row-Standardized Weights
  cat("Attempting Row-Standardized (style='W') conversion...\n")
  
  # Use nb2listw, as it is safer than mat2listw for style="W"
  lw_W <- tryCatch({
    nb2listw(lw_bin$neighbours, style = "W", zero.policy = TRUE)
  }, error = function(e) {
    cat("ERROR in Row-Standardized Conversion:", e$message, "\n")
    return(NULL)
  })
  
  if(!is.null(lw_W)) {
    cat("Row-Standardized conversion successful.\n")
    
    # SANITIZATION CHECK 2
    cat("Checking weights for NaNs...\n")
    has_nan_weights <- any(sapply(lw_W$weights, function(x) any(is.nan(x))))
    if(has_nan_weights) {
      cat("CRITICAL FAILURE: ListW contains NaNs (Division by zero)!\n")
    } else {
      cat("Weights look clean.\n")
      
      # 5. MODEL FIT
      cat("\n--- STEP 5: Model Fitting ---\n")
      # Generate dummy response
      points_df$response <- rnorm(pts_total)
      
      fit <- tryCatch({
        lagsarlm(response ~ intervention, data = points_df, listw = lw_W, zero.policy = TRUE)
      }, error = function(e) {
        cat("ERROR in lagsarlm:", e$message, "\n")
        return(NULL)
      })
      
      if(!is.null(fit)) cat("Model fit successfully!\n")
    }
  }
}

cat("\n--- DIAGNOSTIC COMPLETE ---\n")

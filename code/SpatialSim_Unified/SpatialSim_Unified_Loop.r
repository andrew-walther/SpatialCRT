# ==============================================================================
# Unified Spatial Simulation: Function-Isolated Version
# ==============================================================================

library(tidyverse)
library(sf)
library(spdep)
library(spatialreg)
library(Matrix)

# --- 1. Global Setup ---
set.seed(2024)
ALPHA_TRUE <- 0.2
BETA_TRUE <- 1.0
CELL_SIZE <- 100

# --- 2. The Simulation Engine (One Grid Size at a Time) ---
run_grid_simulation <- function(g_name, rows, cols, n_trt, pts_total) {
  
  cat(sprintf("\n>>> STARTING SIMULATION FOR GRID: %s <<<\n", g_name))
  
  # A. Define Grid Polygons
  bbox <- st_bbox(c(xmin = 0, xmax = cols * CELL_SIZE, ymin = 0, ymax = rows * CELL_SIZE))
  grid_geom <- st_make_grid(st_as_sfc(bbox), n = c(cols, rows), what = "polygons", square = TRUE)
  grid_sf <- st_sf(district = 1:(rows*cols), geometry = grid_geom)
  
  # B. District Adjacency
  nb_dist <- poly2nb(grid_sf, queen = FALSE)
  adj_dist <- nb2mat(nb_dist, style = "B", zero.policy = TRUE)
  
  # C. Combinations & Blocking
  combos_matrix <- combn(1:nrow(grid_sf), n_trt)
  is_block <- apply(combos_matrix, 2, function(x) {
    vec <- rep(0, nrow(grid_sf)); vec[x] <- 1
    int_adj <- any(adj_dist[x, x] > 0)
    ctl_adj <- any(adj_dist[which(vec==0), which(vec==0)] > 0)
    return(!int_adj && !ctl_adj)
  })
  
  # D. Parameters
  sim_params <- expand.grid(
    psi = c(0.5, 0.6, 0.7, 0.8), 
    rho = c(0.00, 0.01), 
    spill = c("TrtNoSpill", "TrtSpill"), 
    stringsAsFactors = FALSE
  )
  
  grid_results <- list()
  
  for (i in 1:nrow(sim_params)) {
    curr_p <- sim_params[i, ]
    cat(sprintf("Condition %d/%d: Psi=%.1f Rho=%.2f %s\n", i, nrow(sim_params), curr_p$psi, curr_p$rho, curr_p$spill))
    
    batch_betas <- numeric(ncol(combos_matrix))
    
    for (c_idx in 1:ncol(combos_matrix)) {
      # 1. Sample Points (Rejection Sampling for Integrity)
      valid_sample <- FALSE
      while(!valid_sample) {
        points_sf <- st_sf(id = 1:pts_total, geometry = st_sample(grid_sf, size = rep(20, nrow(grid_sf))))
        points_sf$district <- as.integer(st_intersects(points_sf, grid_sf))
        if(length(unique(points_sf$district)) == nrow(grid_sf)) valid_sample <- TRUE
      }
      
      # 2. Build Point Weights (Robust Matrix Method)
      dist_to_pts <- split(1:pts_total, points_sf$district)
      W_point_mat <- matrix(0, pts_total, pts_total)
      for (d_idx in 1:nrow(grid_sf)) {
        nbs <- nb_dist[[d_idx]]
        if (all(nbs > 0)) {
          pts_in_d <- dist_to_pts[[as.character(d_idx)]]
          pts_in_nbs <- unlist(dist_to_pts[as.character(nbs)])
          W_point_mat[pts_in_d, pts_in_nbs] <- 1 / length(pts_in_nbs)
        }
      }
      
      # Convert Matrix directly to listw (Safest way to avoid !anyNA error)
      point_listw <- mat2listw(W_point_mat, style = "W", zero.policy = TRUE)
      
      # 3. Treatment & Spillover
      trt_vec <- rep(0, nrow(grid_sf)); trt_vec[combos_matrix[, c_idx]] <- 1
      points_sf$intervention <- trt_vec[points_sf$district]
      point_nb_dist <- nb_dist[points_sf$district]
      has_nb_trt <- sapply(point_nb_dist, function(n) any(trt_vec[n] == 1))
      points_sf$spillover <- if(curr_p$spill == "TrtSpill") as.numeric(has_nb_trt) else as.numeric(has_nb_trt & points_sf$intervention == 0)
      
      # 4. Data Generation
      A_inv <- solve(diag(pts_total) - (curr_p$rho * W_point_mat))
      lin_pred <- ALPHA_TRUE + (BETA_TRUE * points_sf$intervention) + (curr_p$psi * points_sf$spillover) + rnorm(pts_total, 0, 0.1)
      points_sf$response <- as.numeric(A_inv %*% lin_pred)
      
      # 5. Fit Model
      fit <- tryCatch({
        lagsarlm(response ~ intervention + spillover, data = points_sf, listw = point_listw, zero.policy = TRUE)
      }, error = function(e) return(NULL))
      
      batch_betas[c_idx] <- if(!is.null(fit)) coef(fit)[3] else NA
    }
    
    batch_df <- data.frame(Beta = batch_betas, ComboID = 1:ncol(combos_matrix)) %>%
      mutate(Grid = g_name, Psi = curr_p$psi, Rho = curr_p$rho, Spill = curr_p$spill)
    
    grid_results[[i]] <- bind_rows(
      batch_df %>% mutate(Sampling = "SRS"),
      batch_df[is_block, ] %>% mutate(Sampling = "BSS")
    )
  }
  return(bind_rows(grid_results))
}

# --- 3. Run all Geometries ---
res_2x4 <- run_grid_simulation("2x4", 2, 4, 4, 160)
res_3x3 <- run_grid_simulation("3x3", 3, 3, 4, 180)
res_3x4 <- run_grid_simulation("3x4", 3, 4, 6, 240)

# Final Consolidation
final_comprehensive_results <- bind_rows(res_2x4, res_3x3, res_3x4)

##################################################

# --- 4. Statistical Summary Logic ---

# Calculate Metrics for each Condition
analysis_summary <- final_comprehensive_results %>%
  group_by(Grid, Psi, Rho, Spill, Sampling) %>%
  summarise(
    Avg_Estimate = mean(Beta, na.rm = TRUE),
    SD_Estimate  = sd(Beta, na.rm = TRUE),
    Bias         = Avg_Estimate - 1.0, # BETA_TRUE is 1.0
    MSE          = (Bias^2) + (var(Beta, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  # Scale SD by 10 as per your manuscript example
  mutate(SD_Scaled = SD_Estimate * 10)

# Function to create a clean Table 1 style output for a specific grid
create_manuscript_table <- function(grid_name) {
  analysis_summary %>%
    filter(Grid == grid_name) %>%
    # Clean up labels for the table
    mutate(
      Spill_Label = recode(Spill, "TrtNoSpill" = "Control Only", "TrtSpill" = "Intervention & Control"),
      Rho_Label   = sprintf("%s (%.2f)", ifelse(Rho > 0, "Yes", "No"), Rho)
    ) %>%
    select(Spill_Label, Rho_Label, Psi, Sampling, MSE, SD_Scaled, Bias) %>%
    pivot_wider(
      names_from = Sampling, 
      values_from = c(MSE, SD_Scaled, Bias),
      names_glue = "{Sampling}_{.value}"
    ) %>%
    # Reorder columns to match your screenshot: MSE, SD, Bias
    select(Spill_Label, Rho_Label, Psi, 
           SRS_MSE, SRS_SD_Scaled, SRS_Bias,
           BSS_MSE, BSS_SD_Scaled, BSS_Bias) %>%
    arrange(desc(Spill_Label), desc(Rho_Label), Psi)
}

# Generate the three tables
table_2x4 <- create_manuscript_table("2x4")
table_3x3 <- create_manuscript_table("3x3")
table_3x4 <- create_manuscript_table("3x4")

# View the 2x4 table (matches your screenshot)
print("Table 1: 2x4 Spatial Grid Setting")
print(table_2x4)

##################################################

# 1. Boxplot of Estimates across all Grids
ggplot(final_comprehensive_results, aes(x = Sampling, y = Beta, fill = Sampling)) +
  geom_boxplot(outlier.alpha = 0.2) +
  facet_wrap(~Grid) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "Distribution of Estimates by Grid and Design",
       subtitle = "Red line indicates True Beta (1.0)")

# 2. Comparison of Mean MSE (Faceted by Grid)
ggplot(analysis_summary, aes(x = as.factor(Psi), y = MSE, color = Sampling, group = Sampling)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  facet_grid(Grid ~ Spill) +
  theme_bw() +
  labs(title = "MSE Comparison Across Geometries",
       x = "Spillover Effect (Psi)",
       y = "Mean Squared Error (MSE)")

##################################################

# Define metric calculation logic
summary_stats <- final_comprehensive_results %>%
  group_by(Grid, Spill, Rho, Psi, Sampling) %>%
  summarise(
    Mean_Beta = mean(Beta, na.rm = TRUE),
    SD_Raw    = sd(Beta, na.rm = TRUE),
    # MSE = Variance + Bias^2
    MSE       = (sd(Beta, na.rm=TRUE)^2) + (mean(Beta, na.rm=TRUE) - 1.0)^2,
    Bias      = mean(Beta, na.rm = TRUE) - 1.0,
    .groups = "drop"
  ) %>%
  # Apply scaling used in manuscript headers
  mutate(
    SD_Scaled = case_when(
      Grid == "2x4" ~ SD_Raw * 10,  # Manuscript Table 1 [cite: 822]
      Grid == "3x3" ~ SD_Raw * 10,  # Manuscript Table 2 [cite: 890]
      Grid == "3x4" ~ SD_Raw * 100, # Manuscript Table 3 [cite: 955]
      TRUE ~ SD_Raw
    )
  )

# Function to generate a wide Table for a specific grid size
generate_manuscript_table <- function(grid_size) {
  summary_stats %>%
    filter(Grid == grid_size) %>%
    mutate(
      Spill_Label = recode(Spill, "TrtNoSpill" = "Control Only", "TrtSpill" = "Intervention & Control"),
      Rho_Label = ifelse(Rho == 0.01, "Yes (0.01)", "No (0)")
    ) %>%
    select(Spill_Label, Rho_Label, Psi, Sampling, MSE, SD_Scaled, Bias) %>%
    pivot_wider(
      names_from = Sampling,
      values_from = c(MSE, SD_Scaled, Bias),
      names_glue = "{Sampling}_{.value}"
    ) %>%
    # Organize columns to match Tables 1-3 [cite: 822, 890, 955]
    select(Spill_Label, Rho_Label, Psi, 
           SRS_MSE, SRS_SD_Scaled, SRS_Bias, 
           BSS_MSE, BSS_SD_Scaled, BSS_Bias) %>%
    arrange(desc(Spill_Label), desc(Rho_Label), Psi)
}

# View Table 1 (2x4 Grid)
table_1_replicated <- generate_manuscript_table("2x4")
print(table_1_replicated)

# Replicate Figure Plotting (Generalized for any Grid)
plot_manuscript_figure <- function(grid_size, title_num) {
  
  # Filter data for specific grid
  plot_data <- final_comprehensive_results %>% 
    filter(Grid == grid_size) %>%
    mutate(
      Sampling_Label = recode(Sampling, "SRS" = "Random", "BSS" = "Block"),
      Spill_Label = recode(Spill, "TrtNoSpill" = "Control Only", "TrtSpill" = "Intervention & Control"),
      Rho_Expr = factor(Rho, levels = c(0, 0.01), labels = c("rho==0", "rho==0.01")),
      Psi_Fact = as.factor(Psi)
    )
  
  ggplot(plot_data, aes(x = Sampling_Label, y = Beta, fill = Psi_Fact)) +
    geom_boxplot(outlier.size = 0.5) +
    # Create the 2x2 facet layout [cite: 766, 767, 784, 785]
    facet_grid(Spill_Label ~ Rho_Expr, scales = "free_y", labeller = label_parsed) +
    theme_minimal() +
    theme(
      strip.background = element_rect(fill = "white"),
      panel.border = element_rect(color = "grey", fill = NA),
      legend.position = "right"
    ) +
    labs(
      title = paste("Figure", title_num, ": Intervention effect MSE -", grid_size, "grid"),
      x = "Sampling Group",
      y = "Mean Squared Error (MSE)",
      fill = expression(Psi)
    )
}

# Replicate Figures 5, 6, and 7 [cite: 613]
fig_5 <- plot_manuscript_figure("2x4", 5)
fig_6 <- plot_manuscript_figure("3x3", 6)
fig_7 <- plot_manuscript_figure("3x4", 7)

# Display Figure 5
print(fig_5)
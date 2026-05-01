# Load required libraries
library(sf)
library(dplyr)
library(ggplot2)
library(tigris)
library(ggpattern)

# Source the raw data mapping script
source("~/GithubProjects/SpatialCRT/projects/IncidenceDesign/application/cc_mapping_data.R")

#' Fetch and Prepare NC Community College Spatial Data
#'
#' @description Fetches NC shapefiles, projects them, and merges them with 
#' the dynamically sourced primary community college mapping dataset.
#'
#' @param cb Logical. If TRUE, uses lower-resolution boundary files.
#' @return An sf dataframe containing the mapped county data.
#' @export
get_nc_cc_data <- function(cb = TRUE) {
  
  # Fetch NC county shapefiles and transform CRS to NC State Plane (EPSG: 32119)
  nc_counties <- tigris::counties(state = "NC", cb = cb, class = "sf", progress_bar = FALSE) %>%
    sf::st_transform(32119)
  
  # Retrieve the data frame from the sourced script
  cc_mapping <- get_cc_mapping_data()
  
  # Create sequential IDs and Legend Labels for the 58 colleges
  unique_colleges <- sort(unique(cc_mapping$Primary_College))
  college_ids <- data.frame(Primary_College = unique_colleges, ID = 1:length(unique_colleges))
  college_ids$Legend_Label <- paste0(college_ids$ID, ") ", college_ids$Primary_College)
  
  # Merge datasets together
  cc_mapping <- dplyr::left_join(cc_mapping, college_ids, by = "Primary_College")
  nc_cluster_map <- dplyr::left_join(nc_counties, cc_mapping, by = "NAME")
  
  # Establish factor for plotting the shared coverage pattern
  nc_cluster_map$Coverage_Status <- factor(nc_cluster_map$Is_Shared, 
                                           levels = c(FALSE, TRUE), 
                                           labels = c("Exclusive", "Secondary: Roanoke-Chowan CC"))
  return(nc_cluster_map)
}

#' Generate Spatial Map Retaining County Borders
#'
#' @param nc_cluster_map Prepared sf dataframe
#' @return ggplot object
#' @export
generate_county_borders_map <- function(nc_cluster_map) {
  merged_clusters <- nc_cluster_map %>%
    dplyr::group_by(Legend_Label, ID) %>%
    dplyr::summarize(geometry = sf::st_union(geometry), .groups = 'drop')
  
  suppressWarnings({ cluster_centroids <- sf::st_point_on_surface(merged_clusters) })
  
  map_obj <- ggplot2::ggplot() +
    ggpattern::geom_sf_pattern(
      data = nc_cluster_map,
      ggplot2::aes(fill = Legend_Label, pattern = Coverage_Status),
      color = "white", linewidth = 0.2, pattern_fill = "black",   
      pattern_color = NA, pattern_density = 0.1, pattern_spacing = 0.02, pattern_angle = 45        
    ) +
    ggplot2::geom_sf_label(
      data = cluster_centroids, ggplot2::aes(label = ID),
      size = 3.5, fontface = "bold", color = "black",
      fill = ggplot2::alpha("white", 0.75), label.size = 0, label.padding = ggplot2::unit(0.15, "lines")
    ) +
    ggpattern::scale_pattern_manual(
      name = "Overlapping Service", values = c("Exclusive" = "none", "Secondary: Roanoke-Chowan CC" = "stripe")
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        title = "Primary Serving College", ncol = 5, title.position = "top",
        label.theme = ggplot2::element_text(size = 7), keywidth = 0.5, keyheight = 0.5,
        override.aes = list(pattern = "none") 
      )
    ) +
    ggplot2::theme_void() +
    ggplot2::labs(
      title = "North Carolina Community College Service Areas (County Level)",
      subtitle = "Numbered by primary school. Diagonal slashes designate counties co-served by Roanoke-Chowan CC."
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      legend.position = "bottom", legend.box = "vertical",
      legend.title = ggplot2::element_text(face = "bold", size = 10)
    )
  
  return(map_obj)
}

#' Generate Spatial Map with Dissolved Interior Borders
#'
#' @param nc_cluster_map Prepared sf dataframe
#' @return ggplot object
#' @export
generate_dissolved_borders_map <- function(nc_cluster_map) {
  nc_primary_dissolved <- nc_cluster_map %>%
    dplyr::group_by(Legend_Label, ID) %>%
    dplyr::summarize(geometry = sf::st_union(geometry), .groups = 'drop')
  
  suppressWarnings({ cluster_centroids <- sf::st_point_on_surface(nc_primary_dissolved) })
  
  map_obj <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = nc_primary_dissolved, ggplot2::aes(fill = Legend_Label), color = "white", linewidth = 0.4) +
    ggpattern::geom_sf_pattern(
      data = nc_cluster_map, ggplot2::aes(pattern = Coverage_Status),
      fill = NA, color = NA, pattern_fill = "black", pattern_color = NA,       
      pattern_density = 0.1, pattern_spacing = 0.02, pattern_angle = 45        
    ) +
    ggplot2::geom_sf_label(
      data = cluster_centroids, ggplot2::aes(label = ID),
      size = 3.5, fontface = "bold", color = "black",
      fill = ggplot2::alpha("white", 0.75), label.size = 0, label.padding = ggplot2::unit(0.15, "lines")
    ) +
    ggpattern::scale_pattern_manual(
      name = "Overlapping Service", values = c("Exclusive" = "none", "Secondary: Roanoke-Chowan CC" = "stripe")
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        title = "Primary Serving College", ncol = 5, title.position = "top",
        label.theme = ggplot2::element_text(size = 7), keywidth = 0.5, keyheight = 0.5
      )
    ) +
    ggplot2::theme_void() +
    ggplot2::labs(
      title = "North Carolina Community College Service Areas (Contiguous Clusters)",
      subtitle = "Interior borders dissolved. Diagonal slashes designate secondary coverage by Roanoke-Chowan CC."
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      legend.position = "bottom", legend.box = "vertical",
      legend.title = ggplot2::element_text(face = "bold", size = 10)
    )
  
  return(map_obj)
}

# ---------------------------------------------------------
# Execution Block
# ---------------------------------------------------------
spatial_data <- get_nc_cc_data()

map_counties <- generate_county_borders_map(spatial_data)
print(map_counties)

map_dissolved <- generate_dissolved_borders_map(spatial_data)
print(map_dissolved)
#' This script's own directory, captured at source-time
#'
#' @description Self-contained directory detection so `get_cc_mapping_data()`
#' can find `cc_name_crosswalk.csv` regardless of which script sources this
#' file or what working directory is active. Must be evaluated at the top
#' level of this file (not inside a function called later) so `sys.frame(1)`
#' resolves to the frame `source()` creates for this file.
.cc_mapping_script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    normalizePath(dirname(sub("--file=", "", file_arg[1])), mustWork = TRUE)
  } else {
    source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
    if (!is.null(source_file)) {
      normalizePath(dirname(source_file), mustWork = TRUE)
    } else {
      normalizePath(getwd(), mustWork = TRUE)
    }
  }
})

#' Define NC Community College Mapping Data
#'
#' @description Generates the hardcoded 100-county mapping dataframe linking
#' each North Carolina county to its primary serving community college
#' (abbreviated short name), plus the college's full official name joined in
#' from `cc_name_crosswalk.csv`. Both sources independently restrict every
#' county to exactly one college, including Bertie and Northampton, which
#' this study also treats as served solely by Martin CC and Halifax CC
#' respectively (see Roanoke-Chowan CC's restriction to Hertford County).
#'
#' @return A dataframe containing NAME, Primary_College, Full_College_Name,
#'   and Is_Shared.
#' @export
get_cc_mapping_data <- function() {
  mapping <- data.frame(
    NAME = c(
      "Alamance", "Alexander", "Alleghany", "Anson", "Ashe", "Avery", "Beaufort", 
      "Bertie", "Bladen", "Brunswick", "Buncombe", "Burke", "Cabarrus", "Caldwell", 
      "Camden", "Carteret", "Caswell", "Catawba", "Chatham", "Cherokee", "Chowan", 
      "Clay", "Cleveland", "Columbus", "Craven", "Cumberland", "Currituck", "Dare", 
      "Davidson", "Davie", "Duplin", "Durham", "Edgecombe", "Forsyth", "Franklin", 
      "Gaston", "Gates", "Graham", "Granville", "Greene", "Guilford", "Halifax", 
      "Harnett", "Haywood", "Henderson", "Hertford", "Hoke", "Hyde", "Iredell", 
      "Jackson", "Johnston", "Jones", "Lee", "Lenoir", "Lincoln", "Macon", "Madison", 
      "Martin", "McDowell", "Mecklenburg", "Mitchell", "Montgomery", "Moore", "Nash", 
      "New Hanover", "Northampton", "Onslow", "Orange", "Pamlico", "Pasquotank", 
      "Pender", "Perquimans", "Person", "Pitt", "Polk", "Randolph", "Richmond", 
      "Robeson", "Rockingham", "Rowan", "Rutherford", "Sampson", "Scotland", "Stanly", 
      "Stokes", "Surry", "Swain", "Transylvania", "Tyrrell", "Union", "Vance", "Wake", 
      "Warren", "Washington", "Watauga", "Wayne", "Wilkes", "Wilson", "Yadkin", "Yancey"
    ),
    Primary_College = c(
      "Alamance CC", "Catawba Valley CC", "Wilkes CC", "South Piedmont CC", "Wilkes CC", "Mayland CC", "Beaufort County CC",
      "Martin CC", "Bladen CC", "Brunswick CC", "A-B Tech", "Western Piedmont CC", "Rowan-Cabarrus CC", "Caldwell CC&TI",
      "College of The Albemarle", "Carteret CC", "Piedmont CC", "Catawba Valley CC", "Central Carolina CC", "Tri-County CC", "College of The Albemarle",
      "Tri-County CC", "Cleveland CC", "Southeastern CC", "Craven CC", "Fayetteville Tech", "College of The Albemarle", "College of The Albemarle",
      "Davidson-Davie CC", "Davidson-Davie CC", "James Sprunt CC", "Durham Tech", "Edgecombe CC", "Forsyth Tech", "Vance-Granville CC",
      "Gaston College", "College of The Albemarle", "Tri-County CC", "Vance-Granville CC", "Lenoir CC", "Guilford Tech", "Halifax CC",
      "Central Carolina CC", "Haywood CC", "Blue Ridge CC", "Roanoke-Chowan CC", "Sandhills CC", "Beaufort County CC", "Mitchell CC",
      "Southwestern CC", "Johnston CC", "Lenoir CC", "Central Carolina CC", "Lenoir CC", "Gaston College", "Southwestern CC", "A-B Tech",
      "Martin CC", "McDowell Tech", "Central Piedmont", "Mayland CC", "Montgomery CC", "Sandhills CC", "Nash CC",
      "Cape Fear CC", "Halifax CC", "Coastal Carolina CC", "Durham Tech", "Pamlico CC", "College of The Albemarle",
      "Cape Fear CC", "College of The Albemarle", "Piedmont CC", "Pitt CC", "Isothermal CC", "Randolph CC", "Richmond CC",
      "Robeson CC", "Rockingham CC", "Rowan-Cabarrus CC", "Isothermal CC", "Sampson CC", "Richmond CC", "Stanly CC",
      "Forsyth Tech", "Surry CC", "Southwestern CC", "Blue Ridge CC", "Beaufort County CC", "South Piedmont CC", "Vance-Granville CC", "Wake Tech",
      "Vance-Granville CC", "Beaufort County CC", "Caldwell CC&TI", "Wayne CC", "Wilkes CC", "Wilson CC", "Surry CC", "Mayland CC"
    ),
    Is_Shared = c(
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
    )
  )

  crosswalk_path <- file.path(.cc_mapping_script_dir, "..", "data", "cc_name_crosswalk.csv")
  crosswalk <- utils::read.csv(crosswalk_path, stringsAsFactors = FALSE)

  mapping <- merge(mapping, crosswalk,
                    by.x = "Primary_College", by.y = "Short_Name",
                    all.x = TRUE, sort = FALSE)
  names(mapping)[names(mapping) == "Full_Name"] <- "Full_College_Name"

  if (any(is.na(mapping$Full_College_Name))) {
    missing <- unique(mapping$Primary_College[is.na(mapping$Full_College_Name)])
    stop("No full-name crosswalk entry for: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  mapping[, c("NAME", "Primary_College", "Full_College_Name", "Is_Shared")]
}
#' Define NC Community College Mapping Data
#'
#' @description Generates the hardcoded 100-county mapping dataframe linking 
#' each North Carolina county to its primary serving community college.
#'
#' @return A dataframe containing NAME, Primary_College, and Is_Shared.
#' @export
get_cc_mapping_data <- function() {
  data.frame(
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
}
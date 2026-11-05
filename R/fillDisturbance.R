#' Fill missing disturbance values
#'
#' If disturbance values are missing for years in `merge_df` they will be filled
#' from the closest year with data
#'
#' @inheritParams trajectoriesFromNational
#' @param merge_df data.frame. Contains a Year column which will be merged with
#'   `disturbance` to determine missing years to be filled
#'
#' @returns a data.frame with a row for each year in `merge_df`
#' @export
#'
#' @examples
#' disturb <- data.frame(Year = c(2000:2005, 2008:2010), Anthro = 1:9, Fire_excl_anthro = 2)
#' 
#' year_dat <- data.frame(Year = 2000:2010)
#' 
#' fillDisturbance(disturb, year_dat)
#' 
fillDisturbance <- function(disturbance, merge_df){
  if(is.null(disturbance)){
    return(disturbance)
  }
  testTable(disturbance, c("Year", "Anthro", "Fire_excl_anthro"))
  testTable(merge_df, c("Year"))
  
  if(hasName(merge_df, "PopulationName") && hasName(disturbance, "PopulationName")){
    by_use <- c("Year", "PopulationName")
  } else {
    by_use <- "Year"
  }
  
  disturbance <- merge(merge_df, disturbance, by = by_use, all.x = T)
  if (anyNA(select(disturbance, "Anthro", "Fire_excl_anthro"))) {
    warning(
      "Years ",
      filter(disturbance, if_any(c(Anthro, Fire_excl_anthro), is.na)) %>%
        pull(Year) %>% paste0(collapse = ", "),
      " have missing disturbance data. ",
      "Anthro will be filled from the previous year with data and fire will be filled with 0s"
    )
    
    disturbance <- tidyr::fill(disturbance, Anthro, .direction = "downup",
                               .by = any_of("PopulationName")) %>%
      mutate(
        Fire_excl_anthro = tidyr::replace_na(Fire_excl_anthro, 0),
        Total_dist = Fire_excl_anthro + Anthro
      )
    if(anyNA(select(disturbance, "Anthro", "Fire_excl_anthro"))){
      warning("None of the disturbance data is within the requested year range. Setting disturbance to NULL.",
           call. = FALSE)
      disturbance=NULL
    }
  }
  disturbance
}
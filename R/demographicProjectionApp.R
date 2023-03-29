
#' Title
#'
#' @return
#' @export
#'
#' @examples
#' @import caribouMetrics shiny dplyr shinydashboard
#'
demographicProjectionApp <- function(n = 1000) {
  if(rlang::is_installed("BayesianCaribouDemographicProjection")){
    BayesianCaribouDemographicProjection::demographicProjectionApp()
  } else {
    stop("The package: BayesianCaribouDemographicProjection",
         " is required to run the app. \nPlease install it with:",
         "devtools::install_github('LandSciTech/BayesianCaribouDemographicProjection')")
  }
}



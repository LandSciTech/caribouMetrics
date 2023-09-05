
#' Run the Bayesian caribou demographic projection app
#'
#' This requires the "BayesianCaribouDemographicProjection" to be is installed.
#' 
#' @return launches a shiny app in the browser
#' @export
#' @family demography
#'
#' @examples
#' if(FALSE){
#'  demographicProjectionApp()
#' }
#'
demographicProjectionApp <- function() {
  if(rlang::is_installed("BayesianCaribouDemographicProjection")){
    BayesianCaribouDemographicProjection::demographicProjectionApp()
  } else {
    stop("The package: BayesianCaribouDemographicProjection",
         " is required to run the app. \nPlease install it with:",
         "devtools::install_github('LandSciTech/BayesianCaribouDemographicProjection')")
  }
}



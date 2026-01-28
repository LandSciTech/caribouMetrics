#' Projections of population growth from demographic model summaries.
#'
#' @param replicates 
#' @param Rbar,Sbar Mean and standard deviation of R_bar and S_bar over time. See `estimateBayesianRates()$parList` for expected form.
#' @param Riv,Siv Parameters defining the distribution of interannual variation. See `estimateBayesianRates()$parList` for expected form.
#' @param addl_params TO DO explain population mgmt parameters
#' @param type The distribution of interannual variation varies between "beta" or "bbou" model types. 
#' @param doSummary logical. Default TRUE. If FALSE returns unprocessed outcomes from caribouPopGrowth. 
#'  If TRUE returns summaries and (if returnSamples = T) sample trajectories from prepareTrajectories.
#' @param returnSamples logical. If FALSE returns only summaries. If TRUE
#'   returns example trajectories as well. 
#' @param ... Additional arguments passed to `caribouPopGrowth`
#' @return a data.frame
#' @family demography
#'
#' @export
#' 
#' @examples
#'
#'
trajectoriesFromSummary <- function(replicates, N0, Rbar, Sbar, Riv, Siv,  
                  type = "beta", addl_params = list(), doSummary = T, returnSamples = T,
                  nthin=formals(bboutools::bb_fit_survival)$nthin,...){

  # TO DO: check for expected elements of input args
  #Finish documenting trajectoriesFromSummary
  #Tests for trajectoriesFromSummary
  #Update vignettes to use trajectoriesFromSummary, and add variation over time + population mgmt example to vignettes.
  #Think about reducing compute requirements for disturbance projection updating.
  #Eventually - update app to include disturbance, and remove trajectoriesFromSummaryForApp and associated parTab.
  
  if(type=="beta"){
    results <- ratesFromBetaSummary(N0, Rbar, Sbar, Riv, Siv, addl_params, replicates, nthin)
  }else{
    stop("Handle bbou case")
  }
  
  rr <- trajectoriesFromBayesian(results, cPars = addl_params,returnSamples = returnSamples, 
                                 doSummary = doSummary, ...)
  
  return(rr)
  
}


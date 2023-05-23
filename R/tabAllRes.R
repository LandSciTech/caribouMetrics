#' Tabulate results from caribou Bayesian IPM
#' 
#' Create a summary table for each population parameter modelled by [caribouBayesianIPM()]
#'
#'
#' @param rrSurvMod 
#' @param startYear 
#' @param endYear 
#' @param doSummary 
#'
#' @return
#' @export
#'
#' @examples
tabAllRes <- function(rrSurvMod, startYear, endYear, doSummary = T) {
  # rrSurvMod=rr.surv;startYear= minYr;endYear= maxYr
  # rrSurvMod=result;doSummary=T
  
  allParams <- c(
    "S.annual.KM", "R", "Rfemale", "pop.growth", "fpop.size",
    "meanAFsurv", "meanR", "meanRfemale",
    "medianLambda", "meanLambda"
  )
  allParams <- allParams[is.element(allParams, rrSurvMod$parameters.to.save)]
  
  # check rrSurvMod has same number of years as start to end
  if(length(rrSurvMod$BUGSoutput$mean[[rrSurvMod$parameters.to.save[1]]]) !=
     length(startYear:endYear)){
    stop("The model result has different length than startYear:endYear")
  }
  
  
  allResults <- lapply(allParams, getSumStats, rrSurvMod, startYear, endYear,
                       doSummary = doSummary)
  
  allResults <- do.call(rbind, allResults)
  
  allResults <- allResults[order(allResults$Year), ]
  allResults <- allResults[order(allResults$Parameter), ]
  row.names(allResults) <- 1:length(allResults$Year)
  allResults
}

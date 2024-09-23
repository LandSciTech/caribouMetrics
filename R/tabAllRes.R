#' Tabulate results from caribou Bayesian model
#'
#' Create a summary table for each population parameter modelled by
#' [caribouBayesianPM()]
#'
#'
#' @param rrSurvMod rjags model object. Produced by
#'   [caribouBayesianPM()$result]
#' @param startYear integer. Start year of the simulation
#' @param endYear integer. Start year of the simulation
#' @param doSummary logical. Should the results be summarized by year (TRUE) or
#'   should results of each run be returned individually?
#'
#' @return a data.frame. If `doSummary = TRUE` this will contain the mean,
#'   standard deviation and upper and lower credible intervals for all
#'   parameters and years. There is also a probability that the population is
#'   viable calculated as the proportion of runs where the population growth
#'   rate was greater than 0.99. If `doSummary = FALSE` then the data.frame
#'   contains a row for every simulation run and parameter. 
#' @noRd
tabAllRes <- function(rrSurvMod, startYear, endYear, doSummary = T) {
  # rrSurvMod=rr.surv;startYear= minYr;endYear= maxYr
  # rrSurvMod=result;doSummary=T
  
  allParams <- c(
    "S.annual.KM", "R", "Rfemale", "pop.growth", "fpop.size",
    "geomLambda", "meanLambda"
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

#' Get a set of simulation results from the national demographic model
#'
#' @param covTableObs data frame with Anthro,fire_excl_anthro and Year numeric columns. Each is vector of numbers between 0 and 100
#'   representing the percentage of the landscape covered by anthropogenic
#'   disturbance buffered by 500 m, and the percentage covered by fire that does
#'   not overlap anthropogenic disturbance. 
#' @param forceUpdate logical. If the default inputs are used the result is
#'   cached. Set `forceUpdate` to TRUE to ensure the simulations are re-run.
#' @inheritParams demographicCoefficients
#' @inheritParams caribouPopGrowth
#' @param N0 initial population size
#' @param cPars optional. Parameters for calculating composition survey bias term.
#'
#' @return Output from caribouPopGrowth function.
#' 
#' @family demography
#' @export
#'
#' @examples
#' getSimsNational()
getSimsNational <- function(replicates = 1000, N0 = 1000, 
                            covTableObs =  expand.grid(Anthro=seq(0,100,by=1),fire_excl_anthro=0,Year=NA), 
                            useQuantiles  = NULL,
                            populationGrowthTable  = NULL,
                            cPars=getScenarioDefaults(),
                            interannualVar = eval(formals(caribouPopGrowth)$interannualVar)) {
  # replicates=1000;N0=1000;Anthro=seq(0,100,by=1);fire_excl_anthro=0;
  # useQuantiles =NULL;forceUpdate=F
  covTableObs$Total_dist <- covTableObs$Anthro + covTableObs$fire_excl_anthro
  
  if (is.null(populationGrowthTable )) {
    populationGrowthTable  <- caribouMetrics::popGrowthTableJohnsonECCC
  }
  if (is.null(useQuantiles )) {
    popGrowthPars <- demographicCoefficients(
      replicates,
      populationGrowthTable = populationGrowthTable
    )
    rateSamplesAll <- demographicRates(covTable = covTableObs,
                                       popGrowthPars = popGrowthPars,
                                       returnSample = TRUE, useQuantiles = FALSE)
  } else {
    popGrowthPars <- demographicCoefficients(
      replicates, useQuantiles = useQuantiles,
      populationGrowthTable = populationGrowthTable
    )
    rateSamplesAll <- demographicRates(covTable = covTableObs,
                                       popGrowthPars = popGrowthPars,
                                       returnSample = T)
  }
  
  bc = unique(subset(rateSamplesAll,select=replicate));nr=nrow(bc)
  bc$c = compositionBiasCorrection(q=runif(nr,cPars$qMin,cPars$qMax),w=cPars$cowMult,u=runif(nr,cPars$uMin,cPars$uMax),
                                   z=runif(nr,cPars$zMin,cPars$zMax))
  rateSamplesAll$c = NULL; rateSamplesAll= merge(rateSamplesAll, bc)
  
  #print(paste("getSimsNational",mean(bc$c)))

  pars <- merge(data.frame(N0 = N0), rateSamplesAll)
  pars <- cbind(subset(pars,select=-N0), caribouPopGrowth(pars$N0, R_bar = pars$R_bar,
                                       S_bar = pars$S_bar, numSteps = 1,
                                       K = FALSE, c = pars$c, 
                                       interannualVar=interannualVar, progress = FALSE))
  names(pars)[names(pars)=="replicate"]= "id"
  simBig <- pars
  
  return(simBig)
}
#' Summarize results of Bayesian demographic model in tables
#'
#' TODO: Finish documenting
#'
#' @param caribouBayesDemogMod caribou Bayesian demographic model results
#'   produced by calling [caribouBayesianIPM()]
#' @inheritParams caribouBayesianIPM
#' @param oo observed simulation data, produced by calling
#'   [simulateObservations()]
#' @param simBig National simulation results, produced by calling
#'   [getSimsNational()]
#' @param getKSDists logical. Should Kolmogorov–Smirnov distances be calculated?
#'
#' @return a list of tables:
#' * rr.summary.all: Mean parameter values for each year and standard deviation,
#'   upper and lower credible intervals, as well as disturbance and other
#'   parameter input values
#' * 
#' @export
#'
#' @examples
getOutputTables <- function(caribouBayesDemogMod, startYear, endYear, oo, simBig,
                            getKSDists) {
  # result=out$result;startYear=minYr;endYear=maxYr;survInput=out$survInput;
  # oo=oo;simBig=simBig
  
  result <- caribouBayesDemogMod$result
  survInput <- caribouBayesDemogMod$survInput
  
  # get summary info for plots
  rr.summary <- tabAllRes(result, startYear, endYear)
  
  if (!is.element("surv", names(survInput))) {
    if (sum(survInput$event, na.rm = T) > 0) {
      obsSurv <- getKMSurvivalEstimates(survInput)
    } else {
      obsSurv <- unique(subset(survInput, !is.na(enter), select = c(Year)))
      obsSurv$surv <- NA
      obsSurv$Var1 <- obsSurv$Year
    }
  } else {
    obsSurv <- survInput
  }
  
  obsSurv$Mean <- obsSurv$surv
  obsSurv$Year <- as.numeric(gsub("as.factor(Year)=", "", obsSurv$Var1, fixed = T))
  obsSurv <- subset(obsSurv, Year > 1000)
  
  obsSurv$parameter <- "Adult female survival"
  obsSurv$type <- "observed"
  
  trueSurv <- subset(oo$exData, select = c(Year, survival))
  names(trueSurv) <- c("Year", "Mean")
  trueSurv$parameter <- "Adult female survival"
  trueSurv$type <- "true"
  
  obsRec <- subset(oo$ageRatioOut, select = c(Year, Count, Class))
  obsRec <- tidyr::pivot_wider(obsRec, id_cols = c("Year"), names_from = "Class",
                               values_from = "Count")
  obsRec$Mean <- obsRec$calf / obsRec$cow
  obsRec$parameter <- "Recruitment"
  obsRec$type <- "observed"
  
  trueRec <- subset(oo$exData, select = c(Year, recruitment))
  names(trueRec) <- c("Year", "Mean")
  trueRec$parameter <- "Recruitment"
  trueRec$type <- "true"
  
  obsLam <- subset(oo$exData, select = c(Year, lambda))
  names(obsLam) <- c("Year", "Mean")
  obsLam <- movingAveGrowthRate(obsLam, oo$cs$assessmentYrs)
  obsLam$parameter <- "Population growth rate"
  obsLam$type <- "true"
  
  obsSize <- subset(oo$exData, select = c(Year, N))
  names(obsSize) <- c("Year", "Mean")
  # pop size returned from Bayesian model is at the start of the year, not the end.
  obsSize$Year <- obsSize$Year + 1
  obsSize$parameter <- "Female population size"
  obsSize$type <- "true"
  
  obsAll <- rbind(obsLam, obsSize, subset(obsRec, select = names(obsLam)),
                  trueRec, subset(obsSurv, select = names(obsLam)), trueSurv)
  
  simBigO <- subset(simBig$summary, select = c(Anthro, Mean, lower, upper, Parameter))
  names(simBigO) <- c("Anthro", "Mean", "Lower 95% CRI", "Upper 95% CRI", "parameter")
  
  # combine cs and simDisturbance and add to all output tables, nest params in a list
  dist_params <- merge(oo$simDisturbance, oo$cs)
  
  rr.summary <- merge(rr.summary, dist_params)
  simBigO <- merge(simBigO, dist_params)
  obsAll <- merge(obsAll, dist_params)
  rr.summary.all <- rr.summary
  sim.all <- simBigO
  obs.all <- obsAll
  
  if (getKSDists) {
    # get Kolmogorov smirnov distance between samples at each point
    
    variables <- unique(simBig$summary$parameter)
    anthroPts <- unique(subset(rr.summary, select = c(Year, Anthro)))
    # TO DO: make this step faster
    bmSamples <- tabAllRes(result, startYear, endYear, doSummary = F)
    bmSamples$type <- "local"
    
    simSamples <- merge(anthroPts, simBig$samples)
    simSamples$Anthro <- NULL
    simSamples$type <- "national"
    
    allSamples <- rbind(subset(bmSamples,
                               is.element(Parameter, unique(simSamples$Parameter))),
                        simSamples)
    
    ksDists <- allSamples %>%
      group_by(Year, Parameter) %>%
      group_modify(~ {
        getKSDist(.x$Value, .x$type)
      })
  } else {
    ksDists <- unique(subset(rr.summary, select = c(Year, Parameter)))
    ksDists$KSDistance <- NA
    ksDists$KSpvalue <- NA
  }
  return(list(rr.summary.all = rr.summary.all, sim.all = sim.all,
              obs.all = obs.all, ksDists = ksDists))
}

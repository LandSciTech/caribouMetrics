#' Summarize results of Bayesian demographic model in tables
#'
#' Produces summary tables for Bayesian caribou population model
#' results. 
#'
#' @param caribouBayesDemogMod caribou Bayesian demographic model results
#'   produced by calling [caribouBayesianPM()]
#' @inheritParams caribouBayesianPM
#' @param paramTable data.frame. Optional. Scenario parameters see
#'   [simulateObservations()]
#' @param exData data.frame. Optional. Output of [simulateObservations()] that
#'   records the true population metrics of the population that observations
#'   were simulated from.
#' @param simInitial Initial simulation results, produced by calling
#'   [getSimsInitial()]
#'
#' @return a list of tables:
#' * rr.summary.all: Mean parameter values for each year and standard deviation,
#'   upper and lower credible intervals projected by the Bayesian model, as well
#'   as scenario input parameters.
#' * sim.all: Mean parameter values and upper and lower credible intervals from
#'   the initial model for each year, as well as scenario input parameters.
#' * obs.all: Observed parameter values with column "Type" identifying if it is
#'   the "true" value of the simulated population or the "observed" value
#'   simulated based on the collaring program parameters.
#' 
#' @family demography
#' @export
#'
#' @examples
#' scns <- getScenarioDefaults(projYears = 10, obsYears = 10,
#'                             obsAnthroSlope = 1, projAnthroSlope = 5,
#'                             collarCount = 20, cowMult = 5)
#'
#' simO <- simulateObservations(scns)
#'
#' out <- caribouBayesianPM(survData = simO$simSurvObs, recruitData = simO$simRecruitObs,
#'                           disturbance = simO$simDisturbance,
#'                           niters=4)
#'
#' getOutputTables(out, exData = simO$exData, paramTable = simO$paramTable,
#'                 simInitial = getSimsInitial())
                             
getOutputTables <- function(caribouBayesDemogMod, 
                            startYear = min(caribouBayesDemogMod$inData$disturbanceIn$Year), 
                            endYear = max(caribouBayesDemogMod$inData$disturbanceIn$Year), 
                            paramTable = data.frame(param = "observed"),
                            exData = NULL,
                            simInitial = NULL) {
  # caribouBayesDemogMod = out; startYear = oo$minYr;endYear = oo$maxYr; simInitial = simInitial
  # exData = oo$exData; paramTable = oo$paramTable
  
  result <- caribouBayesDemogMod$result
  survInput <- caribouBayesDemogMod$inData$survDataIn
  simObsList <- caribouBayesDemogMod$inData
  
  # get summary info for plots
  rr.summary <- caribouBayesDemogMod$result$summary

  obsSurv <- survInput %>% group_by(PopulationName,Year)%>% summarize(Mortalities = sum(MortalitiesCertain),StartTotal = max(StartTotal)) 
  obsSurv$Mean <- 1-obsSurv$Mortalities/obsSurv$StartTotal
  obsSurv$Parameter <- "Adult female survival"
  obsSurv$MetricTypeID <- "survival"
  obsSurv$Type <- "observed"

  obsRec <- subset(simObsList$recruitDataIn, 
                select = c("PopulationName","Year", "Cows", "Calves"))
  obsRec$Mean <- obsRec$Calves / obsRec$Cows
  obsRec$Parameter <- "Recruitment"
  obsRec$MetricTypeID <- "recruitment"
  obsRec$Type <- "observed"

  if(!is.null(exData)){
    
    exData <- merge(exData,unique(subset(rr.summary,select=c(MetricTypeID,Parameter))),all.x=T)
    names(exData)[names(exData)=="Amount"] = "Mean"
    exData$Type = "true"
  } 
  obsAll <- rbind(subset(obsRec, select = c("PopulationName","Year", "Mean", "Parameter", "MetricTypeID",
                                            "Type")),
                  subset(obsSurv, select = c("PopulationName","Year", "Mean", "Parameter", "MetricTypeID",
                                             "Type")), 
                  subset(exData, select = c("PopulationName","Year", "Mean", "Parameter", "MetricTypeID",
                                             "Type")))
  
  # combine paramTable and simDisturbance and add to all output tables, nest params in a list
  dist_params <- merge(simObsList$disturbanceIn, paramTable)
  
  if(!is.null(simInitial)){
    if(F&&!all(unique(simObsList$disturbanceIn$Anthro) %in% simInitial$summary$Anthro)){
      message("recalculating initial sims to match anthropogenic distubance scenario")
      
      simInitial <- getSimsInitial(Anthro = unique(simObsList$disturbanceIn$Anthro),cPars=paramTable)
    }
    
    simBigO <- merge(simInitial$summary, dist_params)
  } else {
    simBigO <- NULL
  }
  
  rr.summary <- merge(rr.summary, dist_params)
  obsAll <- merge(obsAll, dist_params)

  return(list(rr.summary.all = rr.summary, sim.all = simBigO,
              obs.all = obsAll))
}

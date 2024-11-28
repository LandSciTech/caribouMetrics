#' Summarize results of Bayesian demographic model in tables
#'
#' Produces summary tables for Bayesian caribou population model
#' results. Optionally calculates Kolmogorov–Smirnov distances between Bayesian
#' model results and initial model results.
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
#' @param getKSDists logical. Should Kolmogorov–Smirnov distances be calculated?
#'
#' @return a list of tables:
#' * rr.summary.all: Mean parameter values for each year and standard deviation,
#'   upper and lower credible intervals projected by the Bayesian model, as well
#'   as scenario input parameters.
#' * sim.all: Mean parameter values and upper and lower credible intervals from
#'   the initial model for each year, as well as scenario input parameters.
#' * obs.all: Observed parameter values with column "type" identifying if it is
#'   the "true" value of the simulated population or the "observed" value
#'   simulated based on the collaring program parameters.
#' * kdDists: Kolmogorov–Smirnov distances and p-values, NA if getKSDists is false
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
#' out <- caribouBayesianPM(survData = simO$simSurvObs, ageRatio = simO$ageRatioOut,
#'                           disturbance = simO$simDisturbance,
#'                           Nchains = 1, Niter = 100, Nburn = 10,
#'                           Nthin = 2)
#'
#' getOutputTables(out, exData = simO$exData, paramTable = simO$paramTable,
#'                 simInitial = getSimsInitial(), getKSDists = FALSE)
                             
getOutputTables <- function(caribouBayesDemogMod, 
                            startYear = min(caribouBayesDemogMod$inData$disturbanceIn$Year), 
                            endYear = max(caribouBayesDemogMod$inData$disturbanceIn$Year), 
                            paramTable = data.frame(param = "observed"),
                            exData = NULL,
                            simInitial = NULL,
                            getKSDists = FALSE) {
  # result=out$result;startYear=minYr;endYear=maxYr;survInput=out$survInput;
  # simObsList=oo;simInitial=simBig
  
  result <- caribouBayesDemogMod$result
  survInput <- caribouBayesDemogMod$inData$survDataIn
  simObsList <- caribouBayesDemogMod$inData
  
  # get summary info for plots
  rr.summary <- tabAllRes(result, startYear, endYear)
  
  if (!is.element("surv", names(survInput))) {
    if (sum(survInput$event, na.rm = T) > 0) {
      obsSurv <- getKMSurvivalEstimates(survInput)

    }else{obsSurv=data.frame()} 
    
    
    if(nrow(obsSurv)==0) {
      if(!is.element("id",names(survInput))){
        survInput$id=NA
        survInput$id[!is.na(survInput$enter)]=1
      }
      obsSurv <- unique(subset(survInput, !is.na(survInput$id),
                               select = c("Year")))
      obsSurv$surv <- NA
    }
  } else {
    obsSurv <- survInput
  }
  
  obsSurv$Mean <- obsSurv$surv
  obsSurv$Year <- as.numeric(obsSurv$Year)
  obsSurv <- subset(obsSurv, obsSurv$Year > 1000)
  
  obsSurv$parameter <- "Adult female survival"
  obsSurv$type <- "observed"
  
  obsRec <- subset(simObsList$ageRatioIn, 
                   select = c("Year", "Count", "Class"))
  obsRec <- tidyr::pivot_wider(obsRec, id_cols = c("Year"), names_from = "Class",
                               values_from = "Count")
  obsRec$Mean <- obsRec$calf / obsRec$cow
  obsRec$parameter <- "Recruitment"
  obsRec$type <- "observed"

  if(!is.null(exData)){
    trueSurv <- subset(exData, select = c("Year", "survival"))
    names(trueSurv) <- c("Year", "Mean")
    trueSurv$parameter <- "Adult female survival"
    trueSurv$type <- "true"
    
    trueRec <- subset(exData, select = c("Year", "recruitment"))
    names(trueRec) <- c("Year", "Mean")
    trueRec$parameter <- "Recruitment"
    trueRec$type <- "true"

    trueX <- subset(exData, select = c("Year", "Rfemale"))
    names(trueX) <- c("Year", "Mean")
    trueX$parameter <- "Adjusted recruitment"
    trueX$type <- "true"
    
    obsLam <- subset(exData, select = c("Year", "lambda"))
    names(obsLam) <- c("Year", "Mean")
    obsLam <- movingAveGrowthRate(obsLam, paramTable$assessmentYrs)
    obsLam$parameter <- "Population growth rate"
    obsLam$type <- "true"
    
    obsSize <- subset(exData, select = c("Year", "N"))
    names(obsSize) <- c("Year", "Mean")
    # pop size returned from Bayesian model is at the start of the year, not the end.
    obsSize$Year <- obsSize$Year + 1
    obsSize$parameter <- "Female population size"
    obsSize$type <- "true"
  } else{
    obsLam <- NULL
    obsSize <- NULL
    trueRec <- NULL
    trueX <- NULL
    trueSurv <- NULL
  }
  obsAll <- rbind(obsLam,
                  obsSize,
                  subset(obsRec, select = c("Year", "Mean", "parameter", "type")),
                  trueRec,
                  trueX,
                  subset(obsSurv, select = c("Year", "Mean", "parameter", "type")), 
                  trueSurv)
  
  # combine paramTable and simDisturbance and add to all output tables, nest params in a list
  dist_params <- merge(simObsList$disturbanceIn, paramTable)
  
  if(!is.null(simInitial)){
    if(!all(unique(simObsList$disturbanceIn$Anthro) %in% simInitial$summary$Anthro)){
      message("recalculating initial sims to match anthropogenic distubance scenario")
      
      simInitial <- getSimsInitial(Anthro = unique(simObsList$disturbanceIn$Anthro),cPars=paramTable)
    }
    
    simBigO <- subset(simInitial$summary, select = c("Anthro", "Mean", "lower",
                                                      "upper", "Parameter"))
    names(simBigO) <- c("Anthro", "Mean", "Lower 95% CRI", "Upper 95% CRI", "parameter")
    
    simBigO <- merge(simBigO, dist_params)
  } else {
    simBigO <- NULL
  }
  
  rr.summary <- merge(rr.summary, dist_params)
  obsAll <- merge(obsAll, dist_params)
  rr.summary.all <- rr.summary
  sim.all <- simBigO
  obs.all <- obsAll
  
  if (getKSDists) {
    if(is.null(simInitial)){
      stop("Cannot calculate KSDists without simInitial.", 
           "Set getKSDists = FALSE or provide simInitial.", call. = FALSE)
    }
    # get Kolmogorov smirnov distance between samples at each point
    
    variables <- unique(simInitial$summary$parameter)
    anthroPts <- unique(subset(rr.summary, select = c("Year", "Anthro")))
    # TO DO: make this step faster
    bmSamples <- tabAllRes(result, startYear, endYear, doSummary = F)
    bmSamples$type <- "local"
    
    simSamples <- merge(anthroPts, simInitial$samples)
    simSamples$Anthro <- NULL
    simSamples$type <- "initial"
    
    allSamples <- rbind(subset(bmSamples,
                               is.element(bmSamples$Parameter, unique(simSamples$Parameter))),
                        simSamples)
    
    ksDists <- allSamples %>%
      group_by(.data$Year, .data$Parameter) %>%
      group_modify(~ {
        getKSDist(.x$Value, .x$type)
      })
  } else {
    ksDists <- unique(subset(rr.summary, select = c("Year", "Parameter")))
    ksDists$KSDistance <- NA
    ksDists$KSpvalue <- NA
  }
  return(list(rr.summary.all = rr.summary.all, sim.all = sim.all,
              obs.all = obs.all, ksDists = ksDists))
}

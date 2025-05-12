#' Run the Bayesian population model for multiple parameter sets
#'
#' Define scenarios in a table and [simulateObservations()], run the
#' [caribouBayesianPM()] model and [getOutputTables()] for each scenario.
#'
#' @param scns data.frame. Parameters for the simulations. See
#'   [getScenarioDefaults()] for details.
#' @param ePars list. Additional parameters passed on to
#'   [simulateObservations()]
#' @inheritParams caribouBayesianPM
#' @inheritParams getOutputTables
#' @param printProgress logical. Should the scenario number and parameters be
#'   printed at each step?
#'
#' @return A list similar to [getOutputTables()] where tables for each scenario
#'   have been appended together. Plus an error log for any scenarios that
#'   failed to run.
#'   
#' @family demography
#' @export
#'
#' @examples
#' scns <- expand.grid(
#'   obsYears =c(10, 20), collarCount = c(30, 300), cowMult = 2, collarInterval = 2,
#'   assessmentYrs = 1, iAnthro = 0,
#'   obsAnthroSlope = 0, projAnthroSlope = 0, sQuantile = 0.9,
#'   rQuantile = 0.7, N0 = 1000
#' )
#' 
#' eParsIn <- list(collarOnTime = 1, collarOffTime = 12, collarNumYears = 3)
#' scResults <- runScnSet(scns, eParsIn, getSimsInitial(),
#'                        niters = 10)# only set to speed up example. Normally keep defaults.


runScnSet <- function(scns, ePars, simInitial,
                      printProgress = FALSE,betaPriors="default",niters=formals(bboutools::bb_fit_survival)$niters,nthin=formals(bboutools::bb_fit_survival)$nthin,...) {
  
  # ePars=eParsIn;simInitial=simBig;printProgress=F;niters = formals(bboutools::bb_fit_survival)$niters)
  scns <- getScenarioDefaults(scns)
  errorLog <- list()
  for (p in 1:nrow(scns)) {
    # p=2
    cs <- scns[p, ]
    if (printProgress) {
      print(paste0(c(p, scns[p, ]), collapse = " "))
    }
    
    if(is.element("lQuantile",names(cs))&&!is.na(cs$lQuantile)){
      trajectories <- subset(simInitial$samples,LambdaPercentile == round(cs$lQuantile*100))
    }else{
      trajectories <- simInitial$samples
    }
    trajectories <- subset(trajectories,Replicate==sample(unique(trajectories$Replicate),1))

    oo <- simulateObservations(trajectories, cs, 
                               cowCounts = ePars$cowCounts,
                               freqStartsByYear = ePars$freqStartsByYear,
                               collarNumYears = ePars$collarNumYears,
                               collarOffTime = ePars$collarOffTime,
                               collarOnTime = ePars$collarOnTime,
                               surv_data = simInitial$surv_data, recruit_data=simInitial$recruit_data)
    #plot(plotSurvivalSeries(oo$simSurvObs))
    
    out <- (caribouBayesianPM(
      survData = oo$simSurvObs, recruitData = oo$simRecruitObs,
      disturbance = oo$simDisturbance,
      betaPriors = betaPriors, startYear = oo$minYr, endYear = oo$maxYr,
      N0 = cs$N0,cPars = cs, niters=niters,nthin=nthin,...))
    if (inherits(out, "try-error")) {
      errorLog[[p]] <- list(cs = cs, error = out)
      saveRDS(list(rr.summary.all = rr.summary.all, sim.all = sim.all,
                   obs.all = obs.all, errorLog = errorLog),
              "results/temp.Rds")
      next
    }

    if (inherits(out$result, "try-error")) {
      errorLog[[p]] <- list(cs = cs, error = out$result)
      saveRDS(list(rr.summary.all = rr.summary.all, sim.all = sim.all,
                   obs.all = obs.all, errorLog = errorLog),
              "results/temp.Rds")
      next
    }

    outTabs <- getOutputTables(caribouBayesDemogMod = out, startYear = minYr,
                               endYear = maxYr, simInitial = simInitial,
                               exData = oo$exData, paramTable = oo$paramTable)
    

    if (p == 1) {
      rr.summary.all <- outTabs$rr.summary.all
      sim.all <- outTabs$sim.all
      obs.all <- outTabs$obs.all
    } else {
      rr.summary.all <- rbind(rr.summary.all, outTabs$rr.summary.all)
      sim.all <- rbind(sim.all, outTabs$sim.all)
      obs.all <- rbind(obs.all, outTabs$obs.all)
    }
  }
  if (length(errorLog) > 0) {
    print(errorLog)
  }
  return(list(rr.summary.all = rr.summary.all, sim.all = sim.all,
              obs.all = obs.all, errorLog = errorLog))
}

#' Run the Bayesian population model for multiple parameter sets
#'
#' Define scenarios in a table and [simulateObservations()], run the
#' [caribouBayesianIPM()] model and [getOutputTables()] for each scenario.
#'
#' @param scns data.frame. Parameters for the simulations. See
#'   [getScenarioDefaults()] for details.
#' @param ePars list. Additional parameters passed on to
#'   [simulateObservations()]
#' @inheritParams caribouBayesianIPM
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
#' scResults <- runScnSet(scns, eParsIn, getSimsNational(), getKSDists = FALSE,
#'                        # only set to speed up example. Normally keep defaults.
#'                        Niter = 10, Nburn = 2)


runScnSet <- function(scns, ePars, simNational, survAnalysisMethod = "KaplanMeier",
                      getKSDists = TRUE, printProgress = FALSE, 
                      Niter = formals(caribouBayesianIPM)$Niter,
                      Nburn = formals(caribouBayesianIPM)$Nburn) {
  # ePars=eParsIn;survAnalysisMethod="Exponential";simNational=simBig;getKSDists=T;printProgress=F;Niter = formals(caribouBayesianIPM)$Niter;Nburn = formals(caribouBayesianIPM)$Nburn
  scns <- getScenarioDefaults(scns)
  errorLog <- list()
  for (p in 1:nrow(scns)) {
    # p=1
    cs <- scns[p, ]
    if (printProgress) {
      print(paste0(c(p, scns[p, ]), collapse = " "))
    }

    oo <- simulateObservations(cs, collarNumYears = ePars$collarNumYears,
                               collarOffTime = ePars$collarOffTime,
                               collarOnTime = ePars$collarOnTime)
    
    betaPriors <- getPriors(cs)
    minYr <- min(oo$exData$Year)
    maxYr <- max(oo$simDisturbance$Year)
    out <- try(caribouBayesianIPM(
      survData = oo$simSurvObs, ageRatio = oo$ageRatioOut,
      disturbance = oo$simDisturbance,
      betaPriors = betaPriors, startYear = minYr, endYear = maxYr,
      N0 = oo$exData$N[1], survAnalysisMethod = survAnalysisMethod,
      adjustR = cs$adjustR, assessmentYrs = cs$assessmentYrs, Niter = Niter, 
      Nburn = Nburn
    ))
    if (inherits(out, "try-error")) {
      errorLog[[p]] <- list(cs = cs, error = out)
      saveRDS(list(rr.summary.all = rr.summary.all, sim.all = sim.all,
                   obs.all = obs.all, ksDists = ksDists, errorLog = errorLog),
              "results/temp.Rds")
      next
    }

    if (inherits(out$result, "try-error")) {
      errorLog[[p]] <- list(cs = cs, error = out$result)
      saveRDS(list(rr.summary.all = rr.summary.all, sim.all = sim.all,
                   obs.all = obs.all, ksDists = ksDists, errorLog = errorLog),
              "results/temp.Rds")
      next
    }

    outTabs <- getOutputTables(caribouBayesDemogMod = out, startYear = minYr,
                               endYear = maxYr, simNational = simNational,
                               exData = oo$exData, paramTable = oo$paramTable,
                               getKSDists = getKSDists)

    if (p == 1) {
      rr.summary.all <- outTabs$rr.summary.all
      sim.all <- outTabs$sim.all
      obs.all <- outTabs$obs.all
      ksDists <- merge(outTabs$ksDists, cs)
    } else {
      rr.summary.all <- rbind(rr.summary.all, outTabs$rr.summary.all)
      sim.all <- rbind(sim.all, outTabs$sim.all)
      obs.all <- rbind(obs.all, outTabs$obs.all)
      ksDists <- rbind(ksDists, merge(outTabs$ksDists, cs))
    }
  }
  if (length(errorLog) > 0) {
    print(errorLog)
  }
  return(list(rr.summary.all = rr.summary.all, sim.all = sim.all,
              obs.all = obs.all, ksDists = ksDists, errorLog = errorLog))
}

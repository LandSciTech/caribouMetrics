runScnSet <- function(scns, ePars, simBig, survAnalysisMethod = "KaplanMeier",
                      getKSDists = T, printProgress = F) {
  # ePars=eParsIn;survAnalysisMethod="KaplanMeier";getKSDists=T;printProgress=F
  scns <- fillDefaults(scns)
  errorLog <- list()
  for (p in 1:nrow(scns)) {
    # p=1
    cs <- scns[p, ]
    if (printProgress) {
      print(paste0(c(p, scns[p, ]), collapse = " "))
    }

    if (is.element("cw", names(cs))) {
      ePars$cowCounts$Count <- cs$cw
    }
    oo <- simulateObservations(cs, cowCounts = ePars$cowCounts,
                               freqStartsByYear = ePars$freqStartsByYear,
                               collarNumYears = ePars$collarNumYears,
                               collarOffTime = ePars$collarOffTime,
                               collarOnTime = ePars$collarOnTime)
    betaPriors <- getPriors(cs)
    minYr <- min(oo$exData$Year)
    maxYr <- max(oo$simDisturbance$Year)
    out <- try(runRMModel(
      survData = oo$simSurvObs, ageRatio.herd = oo$ageRatioOut,
      disturbance = oo$simDisturbance,
      betaPriors = betaPriors, startYear = minYr, endYear = maxYr,
      N0 = cs$N0, survAnalysisMethod = survAnalysisMethod,
      adjustR = cs$adjustR, assessmentYrs = cs$assessmentYrs
    ))
    if (inherits(out, "try-error")) {
      errorLog[[p]] <- list(cs = cs, error = out)
      saveRDS(list(rr.summary.all = rr.summary.all, sim.all = sim.all,
                   obs.all = obs.all, ksDists = ksDists, errorLog = errorLog),
              "temp.Rds")
      next
    }

    if (inherits(out$result, "try-error")) {
      errorLog[[p]] <- list(cs = cs, error = out$result)
      saveRDS(list(rr.summary.all = rr.summary.all, sim.all = sim.all,
                   obs.all = obs.all, ksDists = ksDists, errorLog = errorLog),
              "temp.Rds")
      next
    }

    outTabs <- getOutputTables(result = out$result, startYear = minYr,
                               endYear = maxYr, survInput = out$survInput,
                               oo = oo, simBig = simBig, getKSDists = getKSDists)

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

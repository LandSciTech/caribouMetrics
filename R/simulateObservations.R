#' Simulate survival data
#'
#' Simulate caribou survival data. First a true population trajectory is
#' simulated following the national model and a disturbance scenario. Then
#' realistic observations are simulated from this true population based on a
#' collaring program with the given parameters.
#'
#' @param cs list. Parameters for the simulations. See [fillDefaults()] for
#'   details.
#' @param printPlot logical. print a plot of the true population trajectory?
#' @param cowCounts data.frame. Number of cows counted in aerial surveys each
#'   year. Must have 3 columns "Year", "Count", and "Class" where class is "cow"
#'   in all rows
#' @param freqStartsByYear data.frame. Number of collars deployed in each year.
#'   Must have 2 columns "Year" and "numStarts"
#' @param collarNumYears integer. Number of years until collar falls off
#' @param collarOffTime integer. Month that collars fall off. A number from 1
#'   (January) to 12 (December)
#' @param collarOnTime integer. Month that collars are deployed. A number from 1
#'   (January) to 12 (December)
#' @param distScen data.frame. Disturbance scenario. Must have columns "Year",
#'   "Anthro", and "fire_excl_anthro" containing the year, percentage of the landscape
#'   covered by anthropogenic disturbance buffered by 500 m, and the percentage covered
#'   by fire that does not overlap anthropogenic disturbance. See
#'   [disturbanceMetrics()]. If NULL the disturbance scenario is simulated based
#'   on `cs`
#' @inheritParams demographicCoefficients
#' TODO: remove writeFiles option? I think yes. If not need to ask for location
#' @param writeFiles should simSurvObs and ageRatioOut results be saved to csv
#'   files in the tabs folder
#'
#' @return a list with elements:
#'   list(minYr=minYr,maxYr=maxYr,simDisturbance=simDisturbance,simSurvObs=simSurvObs,ageRatioOut=ageRatioOut,exData=popMetricsWide,cs=cs)
#' @export
#'
#' @examples
simulateObservations <- function(cs, printPlot = F, cowCounts,
                                 freqStartsByYear,
                                 collarNumYears = 4, collarOffTime = 5,
                                 collarOnTime = 8, distScen = NULL,
                                 populationGrowthTable = caribouMetrics::popGrowthTableJohnsonECCC,
                                 survivalModelNumber = "M1",
                                 recruitmentModelNumber = "M4", writeFiles = F) {
  # printPlot=T;cowCounts=ePars$cowCounts;freqStartsByYear=ePars$freqStartsByYear;collarNumYears=ePars$collarNumYears;collarOffTime=ePars$collarOffTime;collarOnTime=ePars$collarOnTime
  # distScen = NULL;popGrowthTable = caribouMetrics::popGrowthTableJohnsonECCC;survivalModelNumber = "M1";recruitmentModelNumber = "M4"
  if (is.character(cowCounts)) {
    cowCounts <- read.csv(cowCounts)
  }
  if (is.character(freqStartsByYear)) {
    freqStartsByYear <- read.csv(freqStartsByYear)
  }

  if(!all(vapply(cs, function(x){length(x)==1}, FUN.VALUE = logical(1)))){
    stop("Each element of cs must have length 1", call. = FALSE)
  }

  testTable(cowCounts, c("Year", "Count", "Class"),
            req_vals = list(Year = cs$iYr:(cs$iYr+cs$P-1)),
            acc_vals = list(Class = "cow"))

  testTable(freqStartsByYear, c("Year", "numStarts"),
            req_vals = list(Year = cs$iYr:(cs$iYr+cs$P-1)))

  # Simulate covariate table
  if (is.null(distScen)) {
    covariates <- simCovariates(cs$iA, cs$iF, cs$P + cs$J, cs$aS, cs$aSf, cs$P + 1)
    simDisturbance <- covariates
    simDisturbance$Year <- cs$iYr + simDisturbance$time - 1

    if (writeFiles) {
      write.csv(simDisturbance, "tabs/simDisturbance.csv") # note default file for UI is always the last scenario run
      write.csv(simDisturbance, paste0("tabs/simDisturbance", cs$label, ".csv"))
    }
  } else {
    simDisturbance <- distScen
    simDisturbance$time <- simDisturbance$Year - cs$iYr + 1
    simDisturbance <- filter(simDisturbance, Year <= (cs$iYr + cs$P - 1 + cs$J) &
                               Year >= cs$iYr)
  }

  # simulate true population trajectory
  popMetrics <- simTrajectory(
    numYears = cs$P + cs$J, covariates = simDisturbance,
    popGrowthTable = populationGrowthTable,
    survivalModelNumber = survivalModelNumber,
    recruitmentModelNumber = recruitmentModelNumber,
    recSlopeMultiplier = cs$rS,
    sefSlopeMultiplier = cs$sS, recQuantile = cs$rQ, sefQuantile = cs$sQ, N0 = cs$N0, adjustR = cs$adjustR
  )

  simDisturbance$time <- NULL
  if (printPlot) {
    # TO DO: save info on true population dynamics, add to projection plots for comparison
    base1 <- ggplot2::ggplot(data = popMetrics, ggplot2::aes(
      x = Timestep, y = Amount, colour = Replicate,
      group = Replicate
    )) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~MetricTypeID, scales = "free") +
      ggplot2::xlab("Time") +
      ggplot2::theme(legend.position = "none")
    print(base1)
  }

  popMetricsWide <- tidyr::pivot_wider(popMetrics, id_cols = c(Replicate, Timestep), names_from = MetricTypeID, values_from = Amount)
  popMetricsWide$Year <- cs$iYr + popMetricsWide$Timestep - 1

  exData <- subset(popMetricsWide, (Timestep <= cs$P))

  # Now apply observation process model to get simulated calf:cow and survival data.
  # Use sample sizes in example input data e.g. Eaker

  # reduce sim data tables to length of observations prior to max year
  minYr <- cs$iYr
  maxYr <- cs$iYr + cs$P + cs$J - 1

  # simulate survival data from survival probability.
  if (is.element("ri", names(cs))) {
    freqStartsByYear <- subset(freqStartsByYear, is.element(Year, unique(exData$Year)))
    renewYrs <- intersect(min(exData$Year) + seq(0, 100) * cs$ri, unique(exData$Year))
    freqStartsByYear$numStarts[!is.element(freqStartsByYear$Year, renewYrs)] <- 0
  } else {
    renewYrs <- unique(freqStartsByYear$Year)
  }

  if (is.element("st", names(cs))) {
    freqStartsByYear$numStarts[is.element(freqStartsByYear$Year, renewYrs)] <- cs$st
    simSurvObs <- simSurvivalData(freqStartsByYear, exData, collarNumYears, collarOffTime, collarOnTime, topUp = T)
  } else {
    simSurvObs <- simSurvivalData(freqStartsByYear, exData, collarNumYears, collarOffTime, collarOnTime)
  }
  # if cmult is provided, set cows as a function of number of surviving cows at month 5
  if (is.element("cmult", names(cs))) {
    survsCalving <- subset(simSurvObs, exit >= 6)

    if (nrow(survsCalving) > 0) {
      cowCounts <- as.data.frame(table(survsCalving$Year))
      names(cowCounts) <- c("Year", "Count")
      cowCounts$Year <- as.numeric(as.character(cowCounts$Year))
      cowCounts$Class <- "cow"
      cowCounts$Count <- cs$cmult * cowCounts$Count
    } else {
      cowCounts$Count <- NA
    }
  }
  # given observed total animals & proportion calfs/cows from simulation - get calf/cow ratio
  ageRatioOut <- simCalfCowRatios(cowCounts, minYr, exData)
  ageRatioOut$HerdCode <- "ALAP" # TO DO: remove option for more than one herd in input files, UI, and all associated code...
  ageRatioOut <- subset(ageRatioOut, select = c(names(cowCounts)))
  if (writeFiles) {
    # TO DO: ensure UI code uses column names rather than column positions, and is not sensitive to rearrangement of columns
    write.csv(ageRatioOut, "tabs/simAgeRatio.csv")
    write.csv(ageRatioOut, paste0("tabs/simAgeRatio", cs$label, ".csv"))

    # TO DO: ensure UI code uses column names rather than column positions, and is not sensitive to rearrangement of columns
    write.csv(simSurvObs, "tabs/simSurvData.csv")
    write.csv(simSurvObs, paste0("tabs/simSurvData", cs$label, ".csv"))
    # TO DO: UI option to easily select among available scenarios.
  }
  return(list(minYr = minYr, maxYr = maxYr, simDisturbance = simDisturbance, simSurvObs = simSurvObs, ageRatioOut = ageRatioOut, exData = popMetricsWide, cs = cs))
}

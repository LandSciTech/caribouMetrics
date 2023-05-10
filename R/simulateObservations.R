#' Simulate survival data
#'
#' Simulate caribou survival data. First a true population trajectory is
#' simulated following the national model and a disturbance scenario. Then
#' realistic observations are simulated from this true population based on a
#' collaring program with the given parameters.
#'
#' @param paramTable list. Parameters for the simulations. See [getScenarioDefaults()] for
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
#'   "Anthro", and "fire_excl_anthro" containing the year, percentage of the
#'   landscape covered by anthropogenic disturbance buffered by 500 m, and the
#'   percentage covered by fire that does not overlap anthropogenic disturbance.
#'   See [disturbanceMetrics()]. If NULL the disturbance scenario is simulated
#'   based on `paramTable`
#' @inheritParams demographicCoefficients
#' @param writeFilesDir characater. If not NULL `simSurvObs` and `ageRatioOut`
#'   results will be saved to csv files in the directory provided
#'
#' @return a list with elements:
#'   * minYr: first year in the simulations,
#'   * maxYr: last year in the simulations,
#'   * simDisturbance: a data frame with columns Anthro, fire_excl_anthro, Total_dist, and  Year,
#'   * simSurvObs: a data frame of survival data with columns id, Year, event, enter, and exit,
#'   * ageRatioOut: a data frame of calf cow counts for each year with columns Year, Count, and Class,
#'   * exData: a tibble of expected population metric based on the national model,
#'   * paramTable: a data frame recording the input parameters for the simulation.
#' @export
#'
#' @examples
#' scns <- getScenarioDefaults(projYears = 10, obsYears = 10)
#' simulateObservations(scns,
#'                      freqStartsByYear = data.frame(Year = 2014:2023,
#'                                                    numStarts = 10),
#'                      cowCounts = data.frame(Year = 2014:2023,
#'                                             Count = 10,
#'                                             Class = "cow"))
simulateObservations <- function(paramTable, cowCounts,
                                 freqStartsByYear,
                                 printPlot = FALSE,
                                 collarNumYears = 4, collarOffTime = 5,
                                 collarOnTime = 8, distScen = NULL,
                                 populationGrowthTable = caribouMetrics::popGrowthTableJohnsonECCC,
                                 survivalModelNumber = "M1",
                                 recruitmentModelNumber = "M4",
                                 writeFilesDir = NULL) {
  # printPlot=T;cowCounts=ePars$cowCounts;freqStartsByYear=ePars$freqStartsByYear;
  # collarNumYears=ePars$collarNumYears;collarOffTime=ePars$collarOffTime;
  # collarOnTime=ePars$collarOnTime
  # distScen = NULL;popGrowthTable = caribouMetrics::popGrowthTableJohnsonECCC;
  # survivalModelNumber = "M1";recruitmentModelNumber = "M4"
  if (is.character(cowCounts)) {
    cowCounts <- read.csv(cowCounts)
  }
  if (is.character(freqStartsByYear)) {
    freqStartsByYear <- read.csv(freqStartsByYear)
  }

  if(!all(vapply(paramTable, function(x){length(x)==1}, FUN.VALUE = logical(1)))){
    stop("Each element of paramTable must have length 1", call. = FALSE)
  }

  testTable(cowCounts, c("Year", "Count", "Class"),
            req_vals = list(Year = paramTable$startYear:(paramTable$startYear+paramTable$obsYears-1)),
            acc_vals = list(Class = "cow"))

  testTable(freqStartsByYear, c("Year", "numStarts"),
            req_vals = list(Year = paramTable$startYear:(paramTable$startYear+paramTable$obsYears-1)))

  # Simulate covariate table
  if (is.null(distScen)) {
    covariates <- simCovariates(paramTable$iAnthro, paramTable$iFire, 
                                paramTable$obsYears + paramTable$projYears, 
                                paramTable$obsAnthroSlope, paramTable$projAnthroSlope, 
                                paramTable$obsYears + 1)
    simDisturbance <- covariates
    simDisturbance$Year <- paramTable$startYear + simDisturbance$time - 1

    if (!is.null(writeFilesDir)) {
      write.csv(simDisturbance,
                file.path(writeFilesDir,
                          paste0("simDisturbance", paramTable$label, ".csv")),
                row.names = FALSE)
    }
  } else {
    simDisturbance <- distScen
    simDisturbance$time <- simDisturbance$Year - paramTable$startYear + 1
    simDisturbance <- filter(simDisturbance, Year <= (paramTable$startYear + paramTable$obsYears - 1 + paramTable$projYears) &
                               Year >= paramTable$startYear)
  }

  # simulate true population trajectory
  suppressMessages(
    popMetrics <- simTrajectory(
      numYears = paramTable$obsYears + paramTable$projYears, covariates = simDisturbance,
      popGrowthTable = populationGrowthTable,
      survivalModelNumber = survivalModelNumber,
      recruitmentModelNumber = recruitmentModelNumber,
      recSlopeMultiplier = paramTable$rSlopeMod,
      sefSlopeMultiplier = paramTable$sSlopeMod, recQuantile = paramTable$rQuantile,
      sefQuantile = paramTable$sQuantile,
      N0 = paramTable$N0, adjustR = paramTable$adjustR 
    )
  )

  simDisturbance$time <- NULL
  if (printPlot) {
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

  popMetricsWide <- tidyr::pivot_wider(popMetrics, id_cols = c(Replicate, Timestep),
                                       names_from = MetricTypeID,
                                       values_from = Amount)
  popMetricsWide$Year <- paramTable$startYear + popMetricsWide$Timestep - 1

  exData <- subset(popMetricsWide, (Timestep <= paramTable$obsYears))

  # Now apply observation process model to get simulated calf:cow and survival data.
  # Use sample sizes in example input data e.g. Eaker

  # reduce sim data tables to length of observations prior to max year
  minYr <- paramTable$startYear
  maxYr <- paramTable$startYear + paramTable$obsYears + paramTable$projYears - 1

  # simulate survival data from survival probability.
  if (is.element("collarInterval", names(paramTable))) {
    freqStartsByYear <- subset(freqStartsByYear,
                               is.element(Year, unique(exData$Year)))
    renewYrs <- intersect(min(exData$Year) + seq(0, 100) * paramTable$collarInterval,
                          unique(exData$Year))
    freqStartsByYear$numStarts[!is.element(freqStartsByYear$Year, renewYrs)] <- 0
  } else {
    renewYrs <- unique(freqStartsByYear$Year)
  }

  if (is.element("collarCount", names(paramTable))) {
    freqStartsByYear$numStarts[is.element(freqStartsByYear$Year, renewYrs)] <- paramTable$collarCount
    simSurvObs <- simSurvivalData(freqStartsByYear, exData, collarNumYears,
                                  collarOffTime, collarOnTime, topUp = T)
  } else {
    simSurvObs <- simSurvivalData(freqStartsByYear, exData, collarNumYears,
                                  collarOffTime, collarOnTime)
  }
  # if cowMult is provided, set cows as a function of number of surviving cows at
  # month 5
  if (is.element("cowMult", names(paramTable))) {
    survsCalving <- subset(simSurvObs, exit >= 6)

    if (nrow(survsCalving) > 0) {
      cowCounts <- as.data.frame(table(survsCalving$Year))
      names(cowCounts) <- c("Year", "Count")
      cowCounts$Year <- as.numeric(as.character(cowCounts$Year))
      cowCounts$Class <- "cow"
      cowCounts$Count <- paramTable$cowMult * cowCounts$Count
    } else {
      cowCounts$Count <- NA
    }
  }
  # given observed total animals & proportion calfs/cows from simulation - get
  # calf/cow ratio
  ageRatioOut <- simCalfCowRatios(cowCounts, minYr, exData)

  ageRatioOut <- subset(ageRatioOut, select = c(names(cowCounts)))
  if (!is.null(writeFilesDir)) {
    write.csv(ageRatioOut,
              file.path(writeFilesDir, paste0("simAgeRatio", paramTable$label, ".csv")),
              row.names = FALSE)
    write.csv(simSurvObs,
              file.path(writeFilesDir, paste0("simSurvData", paramTable$label, ".csv")),
              row.names = FALSE)
  }
  return(list(minYr = minYr, maxYr = maxYr, simDisturbance = simDisturbance,
              simSurvObs = simSurvObs, ageRatioOut = ageRatioOut,
              exData = popMetricsWide, paramTable = paramTable))
}

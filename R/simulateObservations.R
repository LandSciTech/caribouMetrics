#' Simulate survival data
#'
#' Simulate caribou survival data. First a true population trajectory is
#' simulated following the national model and a disturbance scenario. Then
#' realistic observations are simulated from this true population based on a
#' caribou monitoring program with the given parameters. Parameters for the
#' caribou monitoring program, disturbance scenario and the true population
#' trajectory can be specified with `getScenarioDefaults()`.
#' 
#' For a detailed description of the process for simulating data see the
#' [vignette](https://landscitech.github.io/caribouMetrics/articles/BayesianDemographicProjection.html#simulation-of-local-population-dynamics-and-monitoring)
#' (`vignette("BayesianDemographicProjection", package = "caribouMetrics")`).
#'
#' @param paramTable data.frame. Parameters for the simulations. See
#'   [getScenarioDefaults()] for details.
#' @param printPlot logical. print a plot of the true population trajectory?
#' @param cowCounts data.frame. Optional. Number of cows counted in aerial
#'   surveys each year. If NULL, and `paramTable` contains `cowMult` the number
#'   of cows that survive calving based on the collar data is multiplied by
#'   `cowMult` to determine the number of cows counted in aerial surveys. If
#'   `paramTable` does not contain `cowMult` `paramTable$cowCount` is used to
#'   set the number of cows counted in aerial surveys each year. If a data.frame
#'   is provided it must have 3 columns "Year", "Count", and "Class" where class
#'   is "cow" in all rows.
#' @param freqStartsByYear data.frame. Optional. Number of collars deployed in
#'   each year. If NULL `paramTable$collarCount` is used as the target number of
#'   collars and each year that collars are deployed they will be topped up to
#'   this number. If a data.frame is provided it must have 2 columns "Year" and
#'   "numStarts" and the "numStarts" is the absolute number of collars deployed
#'   in that year.
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
#'   * exData: a tibble of expected population metrics based on the national model,
#'   * paramTable: a data frame recording the input parameters for the simulation.
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
simulateObservations <- function(paramTable, cowCounts = NULL,
                                 freqStartsByYear = NULL,
                                 printPlot = FALSE,
                                 collarNumYears = 4, collarOffTime = 5,
                                 collarOnTime = 8, distScen = NULL,
                                 populationGrowthTable = caribouMetrics::popGrowthTableJohnsonECCC,
                                 survivalModelNumber = "M1",
                                 recruitmentModelNumber = "M4",
                                 writeFilesDir = NULL) {
  # paramTable=cs;printPlot=T;cowCounts=NULL;freqStartsByYear=NULL;
  # collarNumYears=ePars$collarNumYears;collarOffTime=ePars$collarOffTime;
  # collarOnTime=ePars$collarOnTime
  # distScen = NULL;populationGrowthTable = caribouMetrics::popGrowthTableJohnsonECCC;
  # survivalModelNumber = "M1";recruitmentModelNumber = "M4";writeFilesDir=NULL
  if (is.character(cowCounts)) {
    cowCounts <- read.csv(cowCounts)
  }
  if (is.character(freqStartsByYear)) {
    freqStartsByYear <- read.csv(freqStartsByYear)
  }

  if(!all(vapply(paramTable, function(x){length(x)==1}, FUN.VALUE = logical(1)))){
    stop("Each element of paramTable must have length 1", call. = FALSE)
  }
  
  # if cowCounts and freqStartsByYear not provided build from cowCount and collarCount
  cowCountsIn <- cowCounts
  freqStartsByYearIn <- freqStartsByYear
  
  if(!is.null(cowCounts)){
    testTable(cowCounts, c("Year", "Count", "Class"),
              req_vals = list(Year = paramTable$startYear:(paramTable$startYear+paramTable$obsYears-1)),
              acc_vals = list(Class = "cow"))
  } else if(!is.null(paramTable$cowCount)){
    cowCounts <- data.frame(Year = paramTable$startYear:
                              (paramTable$startYear+paramTable$obsYears-1),
                            Count = paramTable$cowCount,
                            Class = "cow")
  } else if(is.null(paramTable$cowCount) & is.null(paramTable$cowMult)){
    stop("One of cowCounts or paramTable$cowCount must be provided", 
         call. = FALSE)
  }
  
  if(!is.null(freqStartsByYear)){
    testTable(freqStartsByYear, c("Year", "numStarts"),
              acc_vals = list(Year = paramTable$startYear:(paramTable$startYear+paramTable$obsYears-1)))
  } else if(!is.null(paramTable$collarCount)){
    freqStartsByYear <- data.frame(Year = paramTable$startYear:
                                     (paramTable$startYear+paramTable$obsYears-1),
                                   numStarts = paramTable$collarCount)
  }else {
    stop("One of freqStartsByYear or paramTable$collarCount must be provided",
         call. = FALSE)
  }

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
    simDisturbance <- filter(simDisturbance, .data$Year <= (paramTable$startYear + paramTable$obsYears - 1 + paramTable$projYears) &
                               .data$Year >= paramTable$startYear)
  }

  # simulate true population trajectory
  suppressMessages(
    popMetrics <- simTrajectory(
      numYears = paramTable$obsYears + paramTable$projYears, 
      covariates = simDisturbance,
      popGrowthTable = populationGrowthTable,
      survivalModelNumber = survivalModelNumber,
      recruitmentModelNumber = recruitmentModelNumber,
      recSlopeMultiplier = paramTable$rSlopeMod,
      sefSlopeMultiplier = paramTable$sSlopeMod, recQuantile = paramTable$rQuantile,
      sefQuantile = paramTable$sQuantile,
      N0 = paramTable$N0, adjustR = paramTable$adjustR,
      cowMult = ifelse(is.null(paramTable$cowMult), 1, paramTable$cowMult),
      qMin = paramTable$qMin, qMax = paramTable$qMax, uMin = paramTable$uMin,
      uMax = paramTable$uMax, zMin = paramTable$zMin, zMax = paramTable$zMax
    )
  )

  simDisturbance$time <- NULL
  if (printPlot) {
    base1 <- ggplot2::ggplot(data = popMetrics, ggplot2::aes(
      x = .data[["Timestep"]], y = .data[["Amount"]], colour = .data[["Replicate"]],
      group = "Replicate"
    )) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap("MetricTypeID", scales = "free") +
      ggplot2::xlab("Time") +
      ggplot2::theme(legend.position = "none")
    print(base1)
  }

  popMetricsWide <- tidyr::pivot_wider(popMetrics, id_cols = c("Replicate", "Timestep"),
                                       names_from = "MetricTypeID",
                                       values_from = "Amount")
  popMetricsWide$Year <- paramTable$startYear + popMetricsWide$Timestep - 1

  exData <- subset(popMetricsWide, (popMetricsWide$Timestep <= paramTable$obsYears))

  # Now apply observation process model to get simulated calf:cow and survival data.
  # Use sample sizes in example input data e.g. Eaker

  # reduce sim data tables to length of observations prior to max year
  minYr <- paramTable$startYear
  maxYr <- paramTable$startYear + paramTable$obsYears + paramTable$projYears - 1

  # simulate survival data from survival probability.
  if (is.element("collarInterval", names(paramTable)) & is.null(freqStartsByYearIn)) {
    freqStartsByYear <- subset(freqStartsByYear,
                               is.element(freqStartsByYear$Year, unique(exData$Year)))
    renewYrs <- intersect(min(exData$Year) + seq(0, 100) * paramTable$collarInterval,
                          unique(exData$Year))
    freqStartsByYear$numStarts[!is.element(freqStartsByYear$Year, renewYrs)] <- 0
  } 

  if (is.null(freqStartsByYearIn)) {
    #freqStartsByYear$numStarts=0
    
    simSurvObs <- simSurvivalData(freqStartsByYear, exData, collarNumYears,
                                  collarOffTime, collarOnTime, topUp = T)
  } else {
    simSurvObs <- simSurvivalData(freqStartsByYear, exData, collarNumYears,
                                  collarOffTime, collarOnTime)
  }
  # if cowMult is provided, set cows as a function of number of surviving cows at
  # month 5
  if (is.element("cowMult", names(paramTable)) & is.null(cowCountsIn)) {
    survsCalving <- subset(simSurvObs, simSurvObs$exit >= 6)

    if (nrow(survsCalving) > 0) {
      cowCounts <- as.data.frame(table(survsCalving$Year))
      names(cowCounts) <- c("Year", "Count")
      cowCounts$Year <- as.numeric(as.character(cowCounts$Year))
      cowCounts$Class <- "cow"
      cowCounts$Count <- paramTable$cowMult * cowCounts$Count
    } else {
      cowCounts <- data.frame(Year = unique(simSurvObs$Year), 
                              Count = NA, 
                              Class = "cow")
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

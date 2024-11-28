#' Simulate observations
#'
#' Parameters specify a monitoring program that is applied to simulate observations from the example trajectories. 
#' Parameters for the caribou monitoring program, disturbance scenario and the true population
#' trajectory can be specified with `getScenarioDefaults()`.
#' 
#' For a detailed description of the process for simulating data see the
#' [vignette](https://landscitech.github.io/caribouMetrics/articles/BayesianDemographicProjection.html#simulation-of-local-population-dynamics-and-monitoring)
#' (`vignette("BayesianDemographicProjection", package = "caribouMetrics")`).
#'
#' @param trajectories data.frame.
#' @param paramTable data.frame. Parameters for the simulations. See
#'   [getScenarioDefaults()] for details.
#' @param printPlot logical. print a plot of the true population trajectory?
#' @param cowCounts data.frame. Optional. Number of cows counted in aerial
#'   surveys each year. If NULL, and `paramTable` contains `cowMult` the number
#'   of cows that survive calving based on the collar data is multiplied by
#'   `cowMult` to determine the number of cows counted in aerial surveys. If
#'   `paramTable` does not contain `cowMult` `paramTable$cowCount` is used to
#'   set the number of cows counted in aerial surveys each year. If a data.frame
#'   is provided it must have columns "Year"  and "Cows".
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
#' @inheritParams demographicCoefficients
#' @param writeFilesDir character. If not NULL `simSurvObs` and `ageRatioOut`
#'   results will be saved to csv files in the directory provided
#'
#' @return a list with elements:
#'   * minYr: first year in the simulations,
#'   * maxYr: last year in the simulations,
#'   * simSurvObs: a data frame of survival data in bboutools format,
#'   * ageRatioOut: a data frame of recruitment data in bboutools format,
#'   * exData: a tibble of expected population metrics based on the initial model,
#'   * paramTable: a data frame recording the input parameters for the simulation.
#'
#' @family demography
#' @export
#'
#' @examples
#' scns <- getScenarioDefaults(projYears = 10, obsYears = 10,
#'                             collarCount = 20, cowMult = 5)
#'
#' simO <- simulateObservations(scns)
#' 
simulateObservations <- function(trajectories, paramTable,
         cowCounts = NULL,
         freqStartsByYear = NULL,
         collarNumYears = 4, collarOffTime = 4,
         collarOnTime = 4,
         caribouYearStart = 4,
         writeFilesDir = NULL) {
  #paramTable=cs;cowCounts=NULL;freqStartByYear=NULL; collarNumYears = ePars$collarNumYears
  #collarOffTime = ePars$collarOffTime; collarOnTime = ePars$collarOnTime;caribouYearStart=4; writeFileDir=NULL

  includeTimes = seq((paramTable$startYear+paramTable$preYears):
                       (paramTable$startYear+paramTable$preYears+paramTable$obsYears-1))
  includeYears = unique(trajectories$Year)[includeTimes]
  
  popMetrics = subset(trajectories,is.element(Year,includeYears))
  
  includeYears = unique(popMetrics$Year)
  
  exData <- tidyr::pivot_wider(popMetrics, id_cols = c("Replicate", "Year","Timestep","PopulationName"),
                               names_from = "MetricTypeID",
                               values_from = "Amount")
  
  
  # if cowCounts and freqStartsByYear not provided build from cowCount and collarCount
  cowCountsIn <- cowCounts
  freqStartsByYearIn <- freqStartsByYear
  
  if(!is.null(cowCounts)){
    testTable(cowCounts, c("Year", "Cows"),
              req_vals = list(Year = includeYears))
  } else if(hasName(paramTable, "cowCount")){
    cowCounts <- expand.grid(Year = includeYears,
                             Cows = paramTable$cowCount)
  } else if(!hasName(paramTable, "cowCount") & !hasName(paramTable, "cowMult")){
    stop("One of cowCounts or paramTable$cowMult must be provided",
         call. = FALSE)
  }
  
  if(!is.null(freqStartsByYear)){
    testTable(freqStartsByYear, c("Year","numStarts"),
              acc_vals = list(Year = includeYears))
  } else if(!is.null(paramTable$collarCount)){
    freqStartsByYear <- expand.grid(Year = includeYears,numStarts = paramTable$collarCount)
  }else {
    stop("One of freqStartsByYear or paramTable$collarCount must be provided",
         call. = FALSE)
  }
  
  if(!is.element("PopulationName",names(freqStartsByYear))){
    freqStartsByYear=merge(freqStartsByYear,data.frame(PopulationName=unique(exData$PopulationName)))
  }
  
  # simulate survival data from survival probability.
  if (is.element("collarInterval", names(paramTable)) & is.null(freqStartsByYearIn)) {
    renewYrs <- intersect(min(exData$Year) + seq(0, 100) * paramTable$collarInterval,
                          unique(exData$Year))
    freqStartsByYear$numStarts[!is.element(freqStartsByYear$Year, renewYrs)] <- 0
  }
  
  if (is.null(freqStartsByYearIn)) {
    #freqStartsByYear$numStarts=0
    simSurvObs <- simSurvivalData(freqStartsByYear, exData, collarNumYears, collarOffTime,
                                  collarOnTime, caribouYearStart,topUp = T)
  } else {
    simSurvObs <- simSurvivalData(freqStartsByYear, exData, collarNumYears,
                                  collarOffTime, collarOnTime,caribouYearStart)
  }
  
  #plot(plotSurvivalSeries(subset(simSurvObs,Replicate==simSurvObs$Replicate[1])))
  
  
  # if cowMult is provided, set cows as a function of number of surviving cows at
  # year start month
  if (is.element("cowMult", names(paramTable)) & is.null(cowCountsIn)) {
    
    survsCalving <- subset(simSurvObs, simSurvObs$Month == caribouYearStart)
    
    if (nrow(survsCalving) > 0) {
      cowCounts <- subset(survsCalving, select=c("PopulationName","Replicate","Year","StartTotal"))
      cowCounts$Cows <- paramTable$cowMult * cowCounts$StartTotal
    } else {
      cowCounts <- unique(subset(trajectories,select=c("PopulationName","Replicate","Year")))
      cowCounts$Cows <- 0
    }
    cowCounts$Bulls=0;cowCounts$UnknownAdults=0;cowCounts$Yearlings=0
  }
  
  # given observed total animals & proportion calfs/cows from simulation - get
  # calf/cow ratio
  ageRatioOut <- simCalfCowRatios(cowCounts, exData)
  if (!is.null(writeFilesDir)) {
    write.csv(ageRatioOut,
              file.path(writeFilesDir, paste0("simAgeRatio", paramTable$label, ".csv")),
              row.names = FALSE)
    write.csv(simSurvObs,
              file.path(writeFilesDir, paste0("simSurvData", paramTable$label, ".csv")),
              row.names = FALSE)
  }
  
  return(list(simSurvObs = simSurvObs, ageRatioOut = ageRatioOut,
              exData = trajectories, paramTable = paramTable))
}

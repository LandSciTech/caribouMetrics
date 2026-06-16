#' Default scenario parameters.
#'
#' Returns default parameters for scenarios. Use this function to get a combination of [disturbanceDefaults()], [timeDefaults()], [demographyDefaults()], [nationalTrajectoryDefaults()], and [monitoringDefaults()]. If only one of these sets of parameters is needed consider using the relevant component function instead.
#'
#' @param paramTable a data.frame with column names matching the arguments below. Any
#'   columns that are missing will be filled with the default values.
#' @param includeDist logical. Include [disturbanceDefaults()]?
#' @param includeTime logical. Include [timeDefaults()]?
#' @param includeDemography logical. Include [demographyDefaults()]?
#' @param includeNational logical. Include [nationalTrajectoryDefaults()]?
#' @param includeMonitoring logical. Include [monitoringDefaults()]?
#' @inheritParams disturbanceDefaults
#' @inheritParams nationalTrajectoryDefaults
#' @inheritParams monitoringDefaults
#' @inheritParams caribouPopGrowth
#' @inheritParams bayesianTrajectoryWorkflow
#'
#' @return a data.frame of parameter values including a label that combines all
#'   the parameter names and values into a string
#' @export
#' @family demography
#' @examples
#' getScenarioDefaults()
#'
#' # paramTable list takes precedence over argument values
#' getScenarioDefaults(paramTable = data.frame(iFire = 10, iAnthro = 20, obsYears = 1), obsYears = 5)
#' 
getScenarioDefaults <- function(paramTable = NULL,
                         includeDist = T, includeTime = T, includeDemography = T, includeNational = T, includeMonitoring = T,...) {
  
  if(includeDist){
    paramTable <- disturbanceDefaults(paramTable,...)
    
  }else{
    if(includeTime){
      paramTable <- timeDefaults(paramTable,...)
    } 
  }
  
  if(includeNational){
    paramTable <- nationalTrajectoryDefaults(paramTable,...)
  }else{
    if(includeDemography){
      paramTable <- demographyDefaults(paramTable,...)
    }
  }

  if(includeMonitoring){
    paramTable <- monitoringDefaults(paramTable,...)
  }

  # remove columns that are all NA because they should be missing 
  paramTable <- select(paramTable, -where(~all(is.na(.x))))

  paramTable$ID <- seq(1:nrow(paramTable))
  paramTable$label <- ""
  for (n in names(paramTable)[(length(names(paramTable)) - 1):1]) {
    paramTable$label <- paste0(paramTable$label, n, paramTable[[n]], "_")
  }

  return(paramTable)
}

#' Default parameters for specifying scenario durations.
#'
#' Returns default parameter values for scenario durations. 
#' See [simulateObservations()] for additional details.
#'
#' @param paramTable a data.frame with column names matching the arguments below. Any
#'   columns that are missing will be filled with the default values.
#' @param projYears Number of years of projections
#' @param obsYears Number of years of observations
#' @param preYears Number of years before monitoring begins
#' @param curYear year. The current year. All years before are part of the
#'   observation period and years after are part of the projection period.
#' @param startYear year. First year in observation period. Optional, if not provided
#'   it will be calculated from `curYear` and `obsYears`
#'
#' @return a data.frame of parameter values.
#' @export
#' @family demography
#' @examples
#' timeDefaults()
#' 
timeDefaults <- function(paramTable = NULL,
                         projYears=35, obsYears=15, preYears=0, curYear=2023,startYear=NA,...) {
  
  hasYear <- hasName(paramTable,"startYear")|!is.na(startYear)
  paramTable <- addDefaults(paramTable,c(as.list(environment())))
  
  if (!hasYear) {
    paramTable$startYear <- paramTable$curYear - paramTable$obsYears - paramTable$preYears + 1
  }
  
  return(paramTable)
}

#' Default parameters for simulation of disturbance scenarios.
#'
#' Returns default parameter values for disturbance scenarios. 
#' See [simulateObservations()] for additional details.
#'
#' @param paramTable a data.frame with column names matching the arguments below. Any
#'   columns that are missing will be filled with the default values.
#' @param iFire number. Initial fire disturbance percentage.
#' @param iAnthro number. Initial anthropogenic disturbance percentage
#' @param obsAnthroSlope number. Percent change in anthropogenic disturbance per year in the
#'   observation period
#' @param projAnthroSlope number. Percent change in anthropogenic disturbance per year in
#'   the projection period
#' @inheritParams timeDefaults
#'
#' @return a data.frame of parameter values.
#' @export
#' @family demography
#' @examples
#' disturbanceDefaults()
#'
#' # paramTable list takes precedence over argument values
#' disturbanceDefaults(paramTable = data.frame(iFire = 10, iAnthro = 20, obsYears = 1), obsYears = 5)
#' 
disturbanceDefaults <- function(paramTable = NULL,
                                iFire = 0, iAnthro = 0, obsAnthroSlope = 2, projAnthroSlope = 2,...) {
  paramTable <- timeDefaults(paramTable,...)
  return(addDefaults(paramTable,c(as.list(environment()))))
}

#' Default parameters for simulating demographic trajectories.
#'
#' Returns default parameter values for simulating any type of demographic trajectories. 
#' See [trajectoriesFromNational()], [trajectoriesFromBayesian()], [trajectoriesFromSummary()] or [simulateObservations()] for additional details.
#'
#' @param paramTable a data.frame with column names matching the arguments below. Any
#'   columns that are missing will be filled with the default values.
#' @param N0 number. Initial number of adult females.
#' @param qMin number in 0, 1. Minimum ratio of bulls to cows in composition survey groups.
#' @param qMax number in 0, 1. Maximum ratio of bulls to cows in composition survey groups.
#' @param uMin number in 0, 1. Minimum probability of misidentifying young bulls as adult females and vice versa in composition survey.
#' @param uMax number in 0, 1. Maximum probability of misidentifying young bulls as adult females and vice versa in composition survey.
#' @param zMin number in 0, 1. Minimum probability of missing calves in composition survey.
#' @param zMax number in 0, <1. Maximum probability of missing calves in composition survey.
#' @param cowMult number >= 1. The apparent number of adult females per collared animal in composition survey. Set to NA to use `cowCount`.
#' @param lQuantile number in 0, 1. Lambda quantile
#' @param correlateRates logical. Set TRUE to force correlation between recruitment and survival.
#' @return a data.frame of parameter values.
#' @export
#' @family demography
#' @examples
#' demographyDefaults()
#' 
demographyDefaults <- function(paramTable = NULL,
                               N0 = 1000, qMin = 0, qMax = 0, uMin = 0, uMax = 0, zMin = 0, zMax = 0, cowMult = 6,lQuantile = NA, correlateRates = F,...) {
  return(addDefaults(paramTable,c(as.list(environment()))))
}

#' Default parameters for simulating national demographic trajectories.
#'
#' Returns default parameter values for national demographic trajectories. 
#' See [trajectoriesFromNational()] and [simulateObservations()] for additional details.
#'
#' @param paramTable a data.frame with column names matching the arguments below. Any
#'   columns that are missing will be filled with the default values.
#' @param rSlopeMod number. Disturbance-recruitment slope multiplier
#' @param sSlopeMod number. Disturbance-survival slope multiplier
#' @param sQuantile number in 0,1. Survival quantile.
#' @param rQuantile number in 0,1. Recruitment quantile.
#' @inheritParams demographyDefaults
#' @return a data.frame of parameter values.
#' @export
#' @family demography
#' @examples
#' nationalTrajectoryDefaults()
#' 
nationalTrajectoryDefaults <- function(paramTable = NULL,
                                       sQuantile = NA, rQuantile = NA, rSlopeMod = 1, sSlopeMod = 1, interannualVar = list(eval(formals(caribouPopGrowth)$interannualVar)),...) {
  
  paramTable <- demographyDefaults(paramTable,...)
  return(addDefaults(paramTable,c(as.list(environment()))))
}

#' Default parameters for simulating monitoring.
#'
#' Returns default parameter values for monitoring. 
#' See [simulateObservations()] and [bayesianScenariosWorkflow()] for additional details.
#'
#' @param paramTable a data.frame with column names matching the arguments below. Any
#'   columns that are missing will be filled with the default values.
#' @param collarCount number >= 1. The target number of collars active each year. Set to NA to use `freqStartsPerYear` in `simulateObservations()`
#' @param collarInterval number. Optional. Number of years between collar deployments. If
#'   missing assumed to be every year
#' @param cowCount Optional. Only used in `bayesianScenariosWorkflow()` to set the number of cows per
#'   year in recruitment survey
#' @inheritParams demographyDefaults
#' @return a data.frame of parameter values.
#' @export
#' @family demography
#' @examples
#' monitoringDefaults()
#' 
monitoringDefaults <- function(paramTable = NULL,
                               collarInterval=NA, cowCount=NA, collarCount=NA,...) {
  paramTable <- demographyDefaults(paramTable,...)
  
  paramTable <- addDefaults(paramTable,c(as.list(environment())))
  
  if (hasName(paramTable,"cowMult") && any(paramTable$cowMult != 1) && is.element("cowCount", names(paramTable)) && any(!is.na(paramTable$cowCount))) {
    stop("Specify number of cows per year in recruitment survey (cowCount) or",
         " multiplier of number of collared cows in recruitment survey (cowMult) that is different from 1,",
         " but not both.")
  }
  
  if(hasName(paramTable,"cowCount") && sum(paramTable$collarCount>paramTable$N0,na.rm=T)>0){
    warning("Target number of collars collarCount should not exceed initial population size N0.")
  }
  
  if(hasName(paramTable, "collarCount") && 
     hasName(paramTable, "cowMult") && 
     hasName(paramTable, "N0") && 
     sum(paramTable$collarCount*paramTable$cowMult>paramTable$N0,na.rm=T)>0){
    warning("Set cowMult, collarCount and N0 so the expected number of cows in composition surveys does not exceed initial population size N0.")
  }
  
  return(paramTable)
}

addDefaults <- function(paramTable,defList){
  defList$paramTable <- NULL
  if (is.null(paramTable)) {
    paramTable <- as_tibble(defList)
  } else {
    addCols <- as_tibble(defList[which(!names(defList) %in% names(paramTable))])
    addCols <- select(addCols, -where(~all(is.na(.x))))
    if(ncol(addCols)>0){
      # keep all values in paramTable and add any that are missing using values in defList
      paramTable <- bind_cols(paramTable, addCols)
    }
  }

  # order like defList but keep any extra columns not in defList
  paramTable <- select(paramTable, all_of(intersect(names(defList),names(paramTable))), everything())

  return(paramTable)  
}




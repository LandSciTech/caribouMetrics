#' Generate default scenario parameters
#'
#' This function can be used to generate default parameter values used in
#' [simulateObservations()]
#'
#' @param paramTable a data.frame with column names matching the arguments below. Any
#'   columns that are missing will be filled with the default values.
#' @param iFire number. Initial fire disturbance percentage.
#' @param iAnthro number. Initial anthropogenic disturbance percentage
#' @param obsAnthroSlope number. Percent change in anthropogenic disturbance per year in the
#'   observation period
#' @param projAnthroSlope number. Percent change in anthropogenic disturbance per year in
#'   the projection period
#' @param rSlopeMod number. Disturbance-recruitment slope multiplier
#' @param sSlopeMod number. Disturbance-survival slope multiplier
#' @param rQuantile number in 0, 1. Recruitment quantile
#' @param sQuantile number in 0, 1. Survival quantile
#' @param projYears Number of years of projections
#' @param obsYears Number of years of observations
#' @param qMin number in 0, 1. Minimum ratio of bulls to cows in composition survey groups.
#' @param qMax number in 0, 1. Maximum ratio of bulls to cows in composition survey groups.
#' @param uMin number in 0, 1. Minimum probability of misidentifying young bulls as adult females and vice versa in composition survey.
#' @param uMax number in 0, 1. Maximum probability of misidentifying young bulls as adult females and vice versa in composition survey.
#' @param zMin number in 0, 1. Minimum probability of missing calves in composition survey.
#' @param zMax number in 0, <1. Maximum probability of missing calves in composition survey.
#' @param cowMult number >= 1. The apparent number of adult females per collared animal in composition survey.
#' @inheritParams caribouPopGrowth
#' @inheritParams caribouBayesianIPM
#' @param collarInterval number. Optional. Number of years between collar deployments. If
#'   missing assumed to be every year
#' @param cowCount Optional. Only used in `runScnSet()` to set the number of cows per
#'   year in recruitment survey
#' @param collarCount Optional. If provided target number of deployed collars is overwritten with
#'   this value.
#' @param curYear year. The current year. All years before are part of the
#'   observation period and years after are part of the projection period.
#' @param startYear year. First year in observation period. Optional, if not provided
#'   it will be calculated from `curYear` and `obsYears`
#'
#' @references DeCesare, N.J., Hebblewhite, M., Bradley, M., Smith, K.G.,
#'   Hervieux, D. and Neufeld, L., 2012. Estimating ungulate recruitment and
#'   growth rates using age ratios. The Journal of Wildlife Management, 76(1),
#'   pp.144-153.
#'
#'   Eacker, D.R., Hebblewhite, M., Steenweg, R., Russell, M., Flasko, A. and
#'   Hervieux, D., 2019. Web‐based application for threatened woodland caribou
#'   population modeling. Wildlife Society Bulletin, 43(1), pp.167-177.
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
                         iFire = 0, iAnthro = 0, obsAnthroSlope = 2, projAnthroSlope = 2,
                         rSlopeMod = 1, sSlopeMod = 1,
                         rQuantile = 0.5, sQuantile = 0.5, projYears = 35, obsYears = 15, N0 = 1000,
                         adjustR = TRUE, assessmentYrs = 3,qMin=0,qMax =0.6, 
                         uMin = 0, uMax = 0.2, zMin = 0, zMax = 0.2, cowMult = 6,
                         collarInterval = NA, cowCount = NA, 
                         collarCount = NA, startYear = NA,
                         curYear = 2023) {
  defList <- c(as.list(environment()))
  defList$paramTable <- NULL
  if (is.null(paramTable)) {
    paramTable <- as.data.frame(defList)
  } else {
    # keep all values in paramTable and add any that are missing using values in
    # defList
    paramTable <- cbind(paramTable, defList[which(!names(defList) %in% names(paramTable))])
  }

  # remove columns that are all NA because they should be missing and order like
  # defList but keep any extra columns not in defList
  paramTable <- select(paramTable, all_of(names(defList)), everything(), -where(~all(is.na(.x))))

  if (is.element("cowMult", names(paramTable)) & is.element("cowCount", names(paramTable))) {
    stop("Specify number of cows per year in recruitment survey (cowCount) or",
         " multiplier of number of collared cows in recruitment survey (cowMult),",
         " but not both.")
  }
  paramTable$ID <- seq(1:nrow(paramTable))
  paramTable$label <- ""
  for (n in names(paramTable)[(length(names(paramTable)) - 1):1]) {
    paramTable$label <- paste0(paramTable$label, n, paramTable[[n]], "_")
  }

  if (!is.element("startYear", names(paramTable))) {
    paramTable$startYear <- paramTable$curYear - paramTable$obsYears + 1
  }
  return(paramTable)
}

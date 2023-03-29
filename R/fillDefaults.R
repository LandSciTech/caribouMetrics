#' Generate default scenario parameters
#'
#' This function can be used to generate default parameter values used in
#' [simulateObservations()]
#'
#' @param scns a data.frame with column names matching the arguments below. Any
#'   columns that are missing will be filled with the default values.
#' @param iF number. Initial fire disturbance percentage.
#' @param iA number. Initial anthropogenic disturbance percentage
#' @param aS number. Percent change in anthropogenic disturbance per year in the
#'   observation period
#' @param aSf number. Percent change in anthropogenic disturbance per year in
#'   the projection period
#' @param rS number. Disturbance-recruitment slope multiplier
#' @param sS number. Disturbance-survival slope multiplier
#' @param rQ number in 0, 1. Recruitment quantile
#' @param sQ number in 0, 1. Survival quantile
#' @param J Number of years of projections
#' @param P Number of years of observations
#' @param N0 number. Initial population size
#' @param adjustR logical. Adjust R to account for delayed age at first
#'   reproduction (DeCesare et al. 2012; Eacker et al. 2019)?
#' @param assessmentYrs number. #TODO: 
#' @param ri number. Optional. Number of years between collar deployments. If
#'   missing assumed to be every year
#' @param cmult number. Optional. If provided number of cows is `cmult` \*
#'   number of surviving cows at month 5
#' @param cw Optional. Only used in [runScnSet()] to set the number of cows per
#'   year in recruitment survey
#' @param st Optional. If provided number of starts per year is overwritten with
#'   this value.
#' @param curYear year. The current year. All years before are part of the
#'   observation period and years after are part of the projection period.
#' @param iYr year. First year in observation period. Optional, if not provided
#'   it will be calculated from `curYear` and `P`
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
#'
#' @examples
#' fillDefaults()
#'
#' # scns list takes precedence over argument values
#' fillDefaults(scns = data.frame(iF = 10, iA = 20, P = 1), P = 5)
#' 
fillDefaults <- function(scns = NULL,
                         iF = 0, iA = 0, aS = 0, aSf = 4,
                         rS = 1, sS = 1,
                         rQ = 0.5, sQ = 0.5, J = 20, P = 10, N0 = 1000,
                         adjustR = FALSE, assessmentYrs = 1,
                         ri = NA, cmult = NA, cw = NA, st = NA, iYr = NA,
                         curYear = 2023) {
  defList <- c(as.list(environment()))
  defList$scns <- NULL
  if (is.null(scns)) {
    scns <- as.data.frame(defList)
  } else {
    # keep all values in scns and add any that are missing using values in
    # defList
    scns <- cbind(scns, defList[which(!names(defList) %in% names(scns))])
  }

  # remove columns that are all NA because they should be missing and order like
  # defList
  scns <- select(scns, all_of(names(defList)), -where(~all(is.na(.x))))

  if (is.element("cmult", names(scns)) & is.element("cw", names(scns))) {
    stop("Specify number of cows per year in recruitment survey (cw) or",
         " multiplier of number of collared cows in recruitment survey (cmult),",
         " but not both.")
  }
  scns$ID <- seq(1:nrow(scns))
  scns$label <- ""
  for (n in names(scns)[(length(names(scns)) - 1):1]) {
    scns$label <- paste0(scns$label, n, scns[[n]], "_")
  }

  if (!is.element("iYr", names(scns))) {
    scns$iYr <- scns$curYear - scns$P + 1
  }
  return(scns)
}

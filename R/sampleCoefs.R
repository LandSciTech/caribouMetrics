#' Sample coefficients
#'
#' Sample coefficient values from the normal distribution based on the
#' coefficient and standard error
#'
#' @param coefTable data.table. Table of must have columns "Coefficient" for the
#'   name of the coefficient, "Value" for the value of the coefficient and
#'   "StdErr" for the standard error of coefficients
#' @param replicates number of samples to take.
#'
#' @return A list with elements:
#'    \describe{
#'       \item{"coefSamples"}{Bootstrapped coefficients with \code{replicates} 
#'         rows}
#'       \item{"coefValues"}{Coefficient values taken from 
#'         \code{populationGrowthTable}}
#'     }
#'     
#' @examples 
#' cfs <- getCoefs(popGrowthTableJohnsonECCC, "recruitment", "ECCC", "M3")
#' 
#' sampleCoefs(cfs[[1]], 10)
#' 
#' @export

sampleCoefs <- function(coefTable, replicates){
  # Get bootstrapped coefficients
  allCoefs <- coefTable[["Coefficient"]]
  coefSamples <- lapply(X = allCoefs, function(coef){
    if (0) {
      return(rep(as.numeric(coefTable[coefTable$Coefficient == coef, "Value"]), 
                 times = replicates))
    } else {
      vec <- rnorm(n = replicates, 
                   mean = as.numeric(coefTable[coefTable$Coefficient == coef, "Value"]),
                   sd = as.numeric(coefTable[coefTable$Coefficient == coef, "StdErr"]))
      return(vec)
    }
  })
  names(coefSamples) <- allCoefs
  
  coefMatrix <- do.call(cbind, coefSamples)
  
  coefValues <- data.table::data.table(t(coefTable[, "Value"]))
  names(coefValues) <- coefTable[["Coefficient"]]
  
  modList <- list(coefSamples = coefMatrix,
                  coefValues = coefValues)
  return(modList)
}
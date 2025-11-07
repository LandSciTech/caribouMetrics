#' Sample coefficients
#'
#'
#' @param coefTable data.table. Table must have columns "Coefficient" for the
#'   name of the coefficient, "Value" for the value of the coefficient and
#'   "StdErr" for the standard error of coefficients. Typically created with
#'   [subsetNationalCoefs()]
#'
#' @return For `sampleNationalCoefs` a list with elements:
#'   * "coefSamples": Bootstrapped coefficients with `replicates` rows
#'   * "coefValues": Coefficient values taken from `populationGrowthTable`
#'
#' @examples
#' cfs <- subsetNationalCoefs(popGrowthTableJohnsonECCC, "recruitment", "Johnson", "M3")
#'
#' sampleNationalCoefs(cfs[[1]], 10)
#'
#' @rdname getNationalCoefficients
#' @export

sampleNationalCoefs <- function(coefTable, replicates){
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

  coefStdErrs <- data.table::data.table(t(coefTable[, "StdErr"]))
  names(coefStdErrs) <- coefTable[["Coefficient"]]

  modList <- list(coefSamples = coefMatrix,
                  coefValues = coefValues,
                  coefStdErrs=coefStdErrs)
  return(modList)
}

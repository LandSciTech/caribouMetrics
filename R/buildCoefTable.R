#' Create table of coefficients to use in population models, generated using bootstrapping
#'
#' @param caribouCoefTable
#' @param nBootstrap
#' 
#' @export

buildCoefTable <- function(caribouCoefTable, nBootstrap){
  # Get bootstrapped coefficients
  allCoeffs <- caribouCoefTable[["Coefficient"]]
  coeffTable <- lapply(X = allCoeffs, function(coeff){
    if (0) {
      return(rep(as.numeric(caribouCoefTable[caribouCoefTable$Coefficient == coeff, "Value"]), 
                 times = nBootstrap))
    } else {
      vec <- rnorm(n = nBootstrap, 
                   mean = as.numeric(caribouCoefTable[caribouCoefTable$Coefficient == coeff, "Value"]),
                   sd = as.numeric(caribouCoefTable[caribouCoefTable$Coefficient == coeff, "StdErr"]))
      return(vec)
    }
  })
  names(coeffTable) <- allCoeffs
  
  coeffMatrix <- do.call(cbind, coeffTable)
  
  coeffValues <- data.table(t(caribouCoefTable[, "Value"]))
  names(coeffValues) <- caribouCoefTable[["Coefficient"]]
  
  modList <- list(coeffTable = coeffMatrix,
                  coeffValues = coeffValues)
  return(modList)
}
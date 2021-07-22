#' Create and the correct values in a data table used in for the population growth model. 
#' This is a base level function with no defaults set.
#'
#'
#' @param populationGrowthTable
#' @param resVar
#' @param modVer
#' @param modNum
#' 
#' @description ...
#' 
#' @return ...
#' 
#' @export
#' 

makePopDT <- function(populationGrowthTable,
                      resVar,
                      modelVer,
                      modNum){
  ## Check that both a modVer and modNum parameter have been supplied.
  ## Not required for recruitment models
  if (resVar == "femaleSurvival"){
    if (length(modVer) != length(modNum)){
      modVer <- rep(modVer, times = length(modNum))}
  }
  
  testthat::expect_true(length(modVer) == length(modNum), 
                        label = "Please provide one modNum for modVer.
                        length(modVer) == length(modNum)")
  
  Type <- "National"
  
  if (length(Type) != length(modVer)){
    modType <- rep(Type, times = length(modVer))
  } else {
    modType <- Type
  }
  
  modName <- paste(modVer, modNum,Type, sep = "_")
  
  calcFromCI <- function(ci_upper, ci_lower){
    SD <- (ci_upper-ci_lower)/3.92
    return(SD)
  }
  
  DTs <- lapply(seq_along(modVer), FUN = function(modelIndex){
    DT <- populationGrowthTable[responseVariable %in% resVar &
                                  modelVersion  %in% modVer[modelIndex] &
                                  ModelNumber %in% modNum[modelIndex] &
                                  Type %in% modType[modelIndex],]
    if (any(is.na(DT[["StdErr"]]))){ 
      stdErrCalc <- calcFromCI(ci_lower = DT[["lowerCI"]],
                               ci_upper = DT[["upperCI"]])
      stdErrCalc[!is.na(DT[["StdErr"]])]=DT[["StdErr"]][!is.na(DT[["StdErr"]])]
      DT[, StdErr := stdErrCalc]
    }
    return(DT)
  })
  
  names(DTs) <- modName
  
  return(DTs)
}
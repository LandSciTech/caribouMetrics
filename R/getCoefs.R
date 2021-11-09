# Copyright 2021 Tati Micheletti & Her Majesty the Queen in Right of Canada as represented by the Minister of the Environment
# License GPL-3
#NOTICE: This function has been modified from https://github.com/tati-micheletti/caribouPopGrowthModel/blob/master/R/makeDTforPopGrowth.R

#' Get coefficients for model version
#'
#' Get the coefficients for the chosen model version(s), number(s) and response
#' variable from a table of coefficients for many models. If the table does not
#' contain the standard error it is calculated from the confidence interval.
#'
#' @param populationGrowthTable data.frame. \code{\link{popGrowthTableJohnsonECCC}} is
#'   included in the package and should be used in most cases. A custom table of
#'   model parameters can be provided but it must match the column names of
#'   \code{\link{popGrowthTableJohnsonECCC}}.
#' @param resVar character. Response variable, typically "femaleSurvival" or
#'   "recruitment"
#' @param modVer character vector. Which model version(s) to use. Options are
#'   typically "ECCC" for the model used in the ECCC Report (2011) and "Johnson"
#'   for the model used in Johnson et. al. (2020)
#' @param modNum character vector. Which model number(s) to use see
#'   \code{popGrowthTableJohnsonECCC$ModelNumber} for typical options.
#'
#' @return a named list with one element per model version. The names are
#'   \code{modVer_modNum_Type}. Each element contains a data.frame that
#'   is a subset of \code{populationGrowthTable} for the selected model
#'   
#' @examples 
#' getCoefs(popGrowthTableJohnsonECCC, "femaleSurvival", "ECCC", "M1")
#'
#' @export
#' 

getCoefs <- function(populationGrowthTable,
                     resVar,
                     modVer,
                     modNum){
  
  populationGrowthTable <- filter(populationGrowthTable, 
                                  .data$responseVariable == resVar)
  populationGrowthTable <- data.table::as.data.table(populationGrowthTable)
  
  ## Check that both a modVer and modNum parameter have been supplied.
  ## Not required for recruitment models
  if (resVar == "femaleSurvival"){
    if (length(modVer) != length(modNum)){
      modVer <- rep(modVer, times = length(modNum))}
  }
  
  if(length(modVer) != length(modNum)){
    stop("Please provide one modNum for each modVer. length(modVer) == length(modNum)",
         call. = FALSE)
  } 
  selectedMods <- data.frame(modelVersion = modVer, responseVariable = resVar, 
                             ModelNumber = modNum)
  
  missingMods <- anti_join(selectedMods, populationGrowthTable, 
                             by = c("modelVersion", "ModelNumber", "responseVariable"))
  
  if(nrow(missingMods) > 0){
    stop("Model not available. There is no model: ", paste0(missingMods$modelVersion, ", ",
                                                           missingMods$responseVariable, ", ",
                                                           missingMods$ModelNumber, 
                                                           collapse = "\r\n"),
         call. = FALSE)
  }
  
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
    populationGrowthTable = data.table::as.data.table(populationGrowthTable)
    DT <- populationGrowthTable[populationGrowthTable$responseVariable %in% resVar &
                                  populationGrowthTable$modelVersion  %in% modVer[modelIndex] &
                                  populationGrowthTable$ModelNumber %in% modNum[modelIndex] &
                                  populationGrowthTable$Type %in% modType[modelIndex],]
    if (any(is.na(DT[["StdErr"]]))){ 
      DT = data.table::as.data.table(DT)
      stdErrCalc <- calcFromCI(ci_lower = DT[["lowerCI"]],
                               ci_upper = DT[["upperCI"]])
      stdErrCalc[!is.na(DT[["StdErr"]])] = DT[["StdErr"]][!is.na(DT[["StdErr"]])]
      #DT[, StdErr := stdErrCalc]
      DT$StdErr = stdErrCalc
    }
    return(DT)
  })
  
  names(DTs) <- modName
  
  return(DTs)
}
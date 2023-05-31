# Copyright 2021 Tati Micheletti & Her Majesty the Queen in Right of Canada as represented by the Minister of the Environment
# License GPL-3
#NOTICE: This function has been modified from https://github.com/tati-micheletti/caribouPopGrowthModel/blob/master/R/makeDTforPopGrowth.R

#' Get demographic regression model coefficients for model version
#'
#' @param resVar character. Response variable, typically "femaleSurvival" or
#'   "recruitment"
#' @param modNum character vector. Which model number(s) to use see
#'   [popGrowthTableJohnsonECCC] for typical options.
#'
#' @return For `getCoefs`: a named list with one element per model version. The names are
#'   `modelVersion_modNum_Type`. Each element contains a data.frame that is a subset
#'   of `populationGrowthTable` for the selected model
#'
#' @examples
#' getCoefs(popGrowthTableJohnsonECCC, "femaleSurvival", "Johnson", "M1")
#'
#' @rdname demographicCoefficients
#' @export
#' 

getCoefs <- function(populationGrowthTable,
                     resVar,
                     modelVersion,
                     modNum){
  
  populationGrowthTable <- filter(populationGrowthTable, 
                                  .data$responseVariable == resVar)
  populationGrowthTable <- data.table::as.data.table(populationGrowthTable)
  
  ## Check that both a modelVersion and modNum parameter have been supplied.
  ## Not required for recruitment models
  if (resVar == "femaleSurvival"){
    if (length(modelVersion) != length(modNum)){
      modelVersion <- rep(modelVersion, times = length(modNum))}
  }
  
  if(length(modelVersion) != length(modNum)){
    stop("Please provide one modNum for each modelVersion. length(modelVersion) == length(modNum)",
         call. = FALSE)
  } 
  selectedMods <- data.frame(modelVersion = modelVersion, responseVariable = resVar, 
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
  
  if (length(Type) != length(modelVersion)){
    modType <- rep(Type, times = length(modelVersion))
  } else {
    modType <- Type
  }
  
  modName <- paste(modelVersion, modNum,Type, sep = "_")
  
  calcFromCI <- function(ci_upper, ci_lower){
    SD <- (ci_upper-ci_lower)/3.92
    return(SD)
  }
  
  DTs <- lapply(seq_along(modelVersion), FUN = function(modelIndex){
    populationGrowthTable = data.table::as.data.table(populationGrowthTable)
    DT <- populationGrowthTable[populationGrowthTable$responseVariable %in% resVar &
                                  populationGrowthTable$modelVersion  %in% modelVersion[modelIndex] &
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

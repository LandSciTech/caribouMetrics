#' Wrapper function to generate the population growth model predictions under
#' default conditions that match those needed by most users
#'
#' @param covTable
#' @param popGrowthPars
#' @param ignorePrecision
#' @param returnSample
#' @param useQuantiles
#' @param interannualVar
#' 
#' @description ...
#' 
#' @return ...
#' 
#' @export

calcDemoPredictions <- function(covTable,
                                popGrowthPars,
                                ignorePrecision = T,
                                returnSample = F,
                                useQuantiles = F,
                                interannualVar = F){
  
  if (class(interannualVar) == "list") {
    if (length(setdiff(c("Rec","S"), names(interannualVar))) > 0) {
      stop("Expecting interannualVar to be a named list of precision parameters for Rec and S.")
    }
  }
  else {
    if (is.null(interannualVar) || !interannualVar) {
      interannualVar = list(S = F,Rec = F)
    }
    else{
      interannualVar = list(Rec = 1.62 + 3.44, 
                            S= 13.98 + 2.51)      
    }
  }
  
  pred_S <- generatePopGrowthPredictions(covTable = covTable,
                                         coeffTable = popGrowthPars$coeffTable_Survival[["coeffTable"]],
                                         coeffValues = popGrowthPars$coeffTable_Survival[["coeffValues"]],
                                         modelType =  popGrowthPars$modelVersion,
                                         model = "femaleSurvival",
                                         ignorePrecision = ignorePrecision,
                                         returnSample = returnSample,
                                         useQuantiles = useQuantiles,
                                         interannualVar = interannualVar$S)
  
  pred_R <- generatePopGrowthPredictions(covTable = covTable,
                                         coeffTable = popGrowthPars$coeffTable_Recruitment[["coeffTable"]],
                                         coeffValues = popGrowthPars$coeffTable_Recruitment[["coeffValues"]],
                                         modelType =  popGrowthPars$modelVersion,
                                         model = "recruitment",
                                         ignorePrecision = ignorePrecision,
                                         returnSample = returnSample,
                                         useQuantiles = useQuantiles,
                                         interannualVar = interannualVar$Rec)
  rateSamples = pred_S
  names(rateSamples)[names(rateSamples) == "value"] = "S_bar"
  names(rateSamples)[names(rateSamples) == "average"] = "S_bar"
  names(rateSamples)[names(rateSamples) == "stdErr"] = "S_stdErr"
  names(rateSamples)[names(rateSamples) == "PIlow"] = "S_PIlow"
  names(rateSamples)[names(rateSamples) == "PIhigh"] = "S_PIhigh"
  
  rateSamples = merge(rateSamples,pred_R)
  names(rateSamples)[names(rateSamples) == "value"] = "R_bar"
  names(rateSamples)[names(rateSamples) == "average"] = "R_bar"
  names(rateSamples)[names(rateSamples) == "stdErr"] = "R_stdErr"
  names(rateSamples)[names(rateSamples) == "PIlow"] = "R_PIlow"
  names(rateSamples)[names(rateSamples) == "PIhigh"] = "R_PIhigh"
  
  return(rateSamples)  
}
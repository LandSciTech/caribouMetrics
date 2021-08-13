#' Generate default population growth model predictions
#'
#' A wrapper around \code{\link{generatePopGrowthPredictions}} to generate the 
#' population growth model predictions for both survival and recruitment under
#' default conditions that match those needed by most users. 
#'
#' @param covTable data.frame. A table of covariate values to be used in
#'   predictions. Column names must match the coefficient names in
#'   \code{popGrowthTableJohnson}. Each row is a different scenario.
#' @param popGrowthPars list. The result of \code{calcDemoParams}
#' @param ignorePrecision logical. Should the precision of the model be used if
#'   it is available? When precision is used the variance of each population
#'   around the National mean response is considered in addition to the
#'   uncertainty in the coefficient estimates.
#' @param returnSample logical. If TRUE the returned data.frame has replicates *
#'   scenarios rows. If FALSE the returned data.frame has one row per scenario
#'   and additional columns summarizing the variation among replicates. See
#'   Value for details.
#' @param useQuantiles logical. Should the sampling around the National mean use
#'   quantiles or random sampling?
#' @param interannualVar logical or a list with names "S" and "Rec". Should the
#'   interannual variation in population growth be considered? If TRUE a sample
#'   is taken from the beta distribution with precision 5.06 for the recruitment
#'   model and 16.49 for the survival model. If a list the values in "S" and
#'   "Rec" are used for the precision for survival and recruitment,
#'   respectively.
#'
#' @return A data.frame of predictions. The data.frame includes all columns in
#'   \code{covTable} with additional columns depending on \code{returnSample}.
#'   
#'   If \code{returnSample = FALSE} the number of rows is the same as the 
#'   number of rows in \code{covTable}, additional columns are:
#'   \describe{
#'     \item{"S_bar" and "R_bar"}{The mean estimated values of survival and 
#'       recruitment (calves per cow)}
#'     \item{"S_stdErr" and "R_stdErr"}{Standard error of the estimated values}
#'     \item{"S_PIlow"/"S_PIhigh" and "R_PIlow"/"R_PIhigh"}{95% Prediction 
#'       interval for estimated values}
#'   } 
#'   If \code{returnSample = TRUE} the number of rows is \code{nrow(covTable) *
#'    replicates} additional columns are:
#'   \describe{
#'     \item{"scnID"}{A unique identifier for scenarios provided in
#'        \code{covTable}}
#'     \item{"replicate"}{A replicate identifier, unique within each scenario}
#'     \item{"S_bar" and "R_bar"}{The estimated values of survival and 
#'       recruitment (calves per cow)}
#'   } 
#'
#' @seealso \code{\link{generatePopGrowthPredictions}} for more details
#'
#' @export
calcDemoPredictions <- function(covTable,
                                popGrowthPars,
                                ignorePrecision = TRUE,
                                returnSample = FALSE,
                                useQuantiles = FALSE,
                                interannualVar = FALSE){
  
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
  
  if((length(useQuantiles)==1)&&(!useQuantiles)){
    popGrowthPars$coeffTable_S$quantiles=useQuantiles
    popGrowthPars$coeffTable_R$quantiles=useQuantiles
    
  }else{
    q=getQuantiles(nrow(popGrowthPars$coeffTable_S$coeffTable))
    
    if(is.null(popGrowthPars$coeffTable_S$quantiles)){
      popGrowthPars$coeffTable_S$quantiles=sample(q,replace=F)
    }
    if(is.null(popGrowthPars$coeffTable_R$quantiles)){
      popGrowthPars$coeffTable_R$quantiles=sample(q,replace=F)
    }
    
  }
  
  pred_S <- generatePopGrowthPredictions(covTable = covTable,
                                         coeffTable = popGrowthPars$coeffTable_Survival[["coeffTable"]],
                                         coeffValues = popGrowthPars$coeffTable_Survival[["coeffValues"]],
                                         modelType =  popGrowthPars$modelVersion,                                       
                                         useQuantiles = popGrowthPars$coeffTable_S[["quantiles"]],
                                         model = "femaleSurvival",
                                         ignorePrecision = ignorePrecision,
                                         returnSample=returnSample,
                                         interannualVar=interannualVar$S)
  
  pred_R <- generatePopGrowthPredictions(covTable = covTable,
                                         coeffTable = popGrowthPars$coeffTable_Recruitment[["coeffTable"]],
                                         coeffValues = popGrowthPars$coeffTable_Recruitment[["coeffValues"]],
                                         useQuantiles = popGrowthPars$coeffTable_R[["quantiles"]],
                                         modelType =  popGrowthPars$modelVersion,
                                         model = "recruitment",
                                         ignorePrecision = ignorePrecision,
                                         returnSample=returnSample,
                                         interannualVar=interannualVar$Rec)
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
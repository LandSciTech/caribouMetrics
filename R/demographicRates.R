#' Sample expected survival and recruitment rates
#'
#' A wrapper around \code{\link{sampleRates}} to sample 
#' survival and recruitment rates.
#'
#' @param covTable data.frame. A table of covariate values to be used. Column
#'   names must match the coefficient names in \code{\link{popGrowthTableJohnsonECCC}}. 
#'   Each row is a different scenario.
#' @param popGrowthPars list. Coefficient values returned by
#'  \code{demographicCoefficients}
#' @param ignorePrecision logical. Should the precision of the model be used if
#'   it is available? When precision is used variance among populations around the
#'   National mean responses is considered in addition to the uncertainty about 
#'   the coefficient estimates.
#' @param returnSample logical. If TRUE the returned data.frame has replicates *
#'   scenarios rows. If FALSE the returned data.frame has one row per scenario
#'   and additional columns summarizing the variation among replicates. See
#'   Value for details.
#' @param useQuantiles logical. Only relevant when \code{ignorePrecision = FALSE} and
#'   \code{returnSample = TRUE}. If \code{useQuantiles = TRUE}, each replicate population is
#'   assigned to a quantile of the distribution of variation around the expected
#'   values, and remains in that quantile as covariates change. If
#'   \code{useQuantiles = FALSE}, sampling is done independently for each combination of
#'   scenario and replicate, so the value for a particular replicate population
#'   in one scenario is unrelated to the values for that replicate in other
#'   scenarios. If interested in projecting impacts of changing disturbance on
#'   the trajectories of replicate populations set \code{useQuantiles = TRUE}.
#' @param interannualVar logical or a list with names "S" and "Rec". Should the
#'   interannual variation in population growth be considered? If TRUE a sample
#'   is taken from the beta distribution with precision 5.06 for the recruitment
#'   model and 16.49 for the survival model. If a list the values in "S" and
#'   "Rec" are used for the precision for survival and recruitment,
#'   respectively. TO DO: either remove this option from demographicParameters 
#'   function, or update to align with popGrowthJohnson.
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
#'     \item{"S_bar" and "R_bar"}{The expected values of survival and 
#'       recruitment (calves per cow)}
#'   } 
#'
#' @seealso \code{\link{sampleRates}} for more details
#'
#' @export
demographicRates <- function(covTable,
                                popGrowthPars,
                                ignorePrecision = TRUE,
                                returnSample = FALSE,
                                useQuantiles = FALSE,
                                interannualVar = FALSE){
  #covTable=covTableSim; popGrowthPars = popGrowthParsSmall; ignorePrecision = FALSE; returnSample = TRUE;useQuantiles = TRUE;interannualVar = FALSE
  
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
    popGrowthPars$coefSamples_S$quantiles=useQuantiles
    popGrowthPars$coefSamples_R$quantiles=useQuantiles
    
  }else{
    q=getQuantiles(nrow(popGrowthPars$coefSamples_S$coefSamples))
    
    if(is.null(popGrowthPars$coefSamples_S$quantiles)){
      popGrowthPars$coefSamples_S$quantiles=sample(q,replace=F)
    }
    if(is.null(popGrowthPars$coefSamples_R$quantiles)){
      popGrowthPars$coefSamples_R$quantiles=sample(q,replace=F)
    }
    
  }
  
  pred_S <- sampleRates(covTable = covTable,
                        coefSamples = popGrowthPars$coefSamples_Survival[["coefSamples"]],
                        coefValues = popGrowthPars$coefSamples_Survival[["coefValues"]],
                        modVer =  popGrowthPars$modelVersion,                                       
                        useQuantiles = popGrowthPars$coefSamples_S[["quantiles"]],
                        resVar = "femaleSurvival",
                        ignorePrecision = ignorePrecision,
                        returnSample = returnSample,
                        interannualVar = interannualVar$S)
  
  pred_R <- sampleRates(covTable = covTable,
                        coefSamples = popGrowthPars$coefSamples_Recruitment[["coefSamples"]],
                        coefValues = popGrowthPars$coefSamples_Recruitment[["coefValues"]],
                        useQuantiles = popGrowthPars$coefSamples_R[["quantiles"]],
                        modVer =  popGrowthPars$modelVersion,
                        resVar = "recruitment",
                        ignorePrecision = ignorePrecision,
                        returnSample = returnSample,
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
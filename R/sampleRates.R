#' Sample demographic rates
#'
#' Sample expected survival or recruitment rates based on samples of coefficient
#' values and optionally the model precision.
#'
#' \code{coefSamples} and \code{coefValues} can be created with
#' \code{\link{sampleCoefs}}
#'
#' @param coefSamples matrix. Bootstrapped coefficients with one row per
#'   replicate and one column per coefficient
#' @param coefValues data.table. One row table with expected values for each
#'   coefficient
#' @param modVer character. Which model version is being used. Options are
#'   currently, "ECCC" for the model used in the ECCC Report (2011) and
#'   "Johnson" for the model used in Johnson et. al. (2020)
#' @param resVar character. Response variable, typically "femaleSurvival" or
#'   "recruitment"
#' @param ignorePrecision logical. Should the precision of the model be used if
#'   it is available? When precision is used variation among populations around the
#'   National mean responses is considered in addition to the uncertainty about 
#'   the coefficient estimates.
#' @param returnSample logical. If TRUE the returned data.frame has replicates *
#'   scenarios rows. If FALSE the returned data.frame has one row per scenario
#'   and additional columns summarizing the variation among replicates. See
#'   Value for details.
#' @param quantilesToUse numeric vector of length \code{coefSamples}. Only relevant when
#'   \code{ignorePrecision = FALSE}. Each replicate population is assigned to a quantile 
#'   of the distribution of variation around the expected values, and remains in that 
#'   quantile as covariates change. If NULL sampling is done independently for each combination 
#'   of scenario and replicate, so the value for a particular replicate population in one
#'   scenario is unrelated to the values for that replicate in other scenarios.
#'   Useful for projecting impacts of changing disturbance on the
#'   trajectories of replicate populations.
#' @param predInterval numeric vector with length 2. The default 95\% interval is
#'   (\code{c(0.025,0.975)}). Only relevant when \code{returnSample = TRUE} and \code{quantilesToUse = NULL}. 
#'
#' @inheritParams demographicRates
#'
#' @return A data.frame of predictions. The data.frame includes all columns in
#'   \code{covTable} with additional columns depending on \code{returnSample}.
#'   
#'   If \code{returnSample = FALSE} the number of rows is the same as the 
#'   number of rows in \code{covTable}, additional columns are:
#'   \describe{
#'     \item{"average"}{The mean estimated values of the response variable)}
#'     \item{"stdErr"}{Standard error of the estimated values}
#'     \item{"PIlow"/"PIhigh"}{If \code{quantilesToUse = NULL} these are the percentiles given by predInterval. 
#'       If using quantiles, maximum and minimum values are returned.}
#'   } 
#'   If \code{returnSample = TRUE} the number of rows is \code{nrow(covTable) *
#'    replicates} additional columns are:
#'   \describe{
#'     \item{"scnID"}{A unique identifier for scenarios provided in
#'        \code{covTable}}
#'     \item{"replicate"}{A replicate identifier, unique within each scenario}
#'     \item{value}{The expected values of the response variable}
#'   } 
#'
#' @export

sampleRates <- function(covTable,
                        coefSamples,
                        coefValues,
                        modVer,
                        resVar,
                        ignorePrecision,
                        returnSample,
                        quantilesToUse=NULL,
                        predInterval = c(0.025,0.975)){
  
  tictoc::tic(paste0("Elapsed time for caribou prediction for ",
             resVar, " for ", modVer,":"))
  
  whichCovariates <- names(coefValues)[!names(coefValues) %in% c("Intercept", 
                                                                   "intercept",
                                                                   "precision",
                                                                   "Precision")]
  
  missingCovs <- setdiff(whichCovariates, colnames(covTable))
  
  if(length(missingCovs) > 0){
    stop("Covariates missing in covTable: ", paste0(missingCovs, collapse = ", "), 
         call. = FALSE)
  }
  
  covTableRed <- covTable[, whichCovariates]
  
  if (grepl(x = modVer, pattern = "Johnson")) {
    intt <- coefSamples[, which(colnames(coefSamples) %in% c("Intercept", 
                                                           "intercept"))]
    phi <- coefSamples[, which(colnames(coefSamples) %in% c("Precision", 
                                                          "precision"))]
    
    predictedTableSD <- as.matrix(covTableRed) %*%
      t(coefSamples[,-which(colnames(coefSamples) %in% c("Intercept",
                                                       "intercept",
                                                       "precision",
                                                       "Precision"))])
    
    predictedTableSD <- t(apply(predictedTableSD, 1, function(x,intt){
      exp(intt + x)}, intt = intt))

    #predictedTableSD is bootstrap sample of mean. 
    #Precision parameter gives info about the amount of scatter around the mean. 
    #Easier to understand the principles by thinking about a simpler linear regression example. 
    #Using the notation in the regression analysis section of this wikipedia article, predictedTableSD is sample from distribution of yhat_d, rather than y_d: https://en.wikipedia.org/wiki/Prediction_interval 
    if(!ignorePrecision){
      if(length(phi)==0){
        stop("Missing precision parameter. Set ignorePrecision = TRUE or",
             " add a precision column to coefSamples", call. = FALSE)
      }
      if(nrow(predictedTableSD)<1){
        stop("This code assumes at least one row.")
      }
      
      predictedTableSD = t(apply(predictedTableSD,1,betaSample,phi=phi,quantilesToUse=quantilesToUse))
    }
    
    #Note: more relevant measure of uncertainty is perhaps bootstrap prediction interval
    #https://stats.stackexchange.com/questions/226565/bootstrap-prediction-interval
    #When precision is included points are not Gaussian distributed, and SD isn't an overly useful summary metric. 
    
    # Uncertainty across replicates
    predictedSD <- matrixStats::rowSds(predictedTableSD)
    
    if(is.null(quantilesToUse)){
      qqs <- matrixStats::rowQuantiles(predictedTableSD,probs = predInterval)
    }else{
      qqs <- matrixStats::rowQuantiles(predictedTableSD,probs = c(0,1))
    }  
    # Now the model calculations
    predicted <- as.numeric(exp(as.matrix(coefValues)[
      which(colnames(coefValues) %in% c("Intercept", "intercept"))] +
        as.matrix(covTableRed) %*% 
        as.matrix(coefValues)[-which(colnames(coefValues) %in% 
                                        c("Intercept", 
                                          "intercept",
                                          "precision",
                                          "Precision"))]))
  } 
  else {
    if (grepl(x = modVer, pattern = "ECCC")) {
      
      intt <- coefSamples[, which(colnames(coefSamples) %in% c("Intercept", "intercept"))]
      phi <- coefSamples[, which(colnames(coefSamples) %in% c("Precision", 
                                                            "precision"))]
      
      predictedTableSD <- as.matrix(covTableRed) %*%
        t(coefSamples[,-which(colnames(coefSamples) %in% c("Intercept","intercept","precision","Precision"))])
      
      predictedTableSD<-t(apply(predictedTableSD,1,function(x,intt){(intt+x)},intt=intt))
      
      if(resVar == "recruitment"){
        predictedTableSD=predictedTableSD/100
      }
      
      if(!ignorePrecision){
        if(length(phi)==0){
          stop("Missing precision parameter. Set ignorePrecision = TRUE or",
               " add a precision column to coefSamples", call. = FALSE)
        }
        if(nrow(predictedTableSD)<1){
          stop("This code assumes at least one row.")
        }
        
        pp = apply(predictedTableSD,1,normalSample,sd=phi,quantilesToUse=quantilesToUse)
        
        predictedTableSD=t(pp)
        
      }
      
      # Uncertainty across replicates
      predictedSD <- matrixStats::rowSds(predictedTableSD)
      
      if(is.null(quantilesToUse)){
        qqs <- matrixStats::rowQuantiles(predictedTableSD, probs = predInterval)
      }else{
        qqs <- matrixStats::rowQuantiles(predictedTableSD, probs = c(0,1))
      }  
      
      # Now the model calculations
      predicted <- as.numeric((as.matrix(coefValues)[which(colnames(coefValues) %in% c("Intercept",
                                                                                         "intercept"))] +
                                 as.matrix(covTableRed) %*%
                                 as.matrix(coefValues)[-which(colnames(coefValues) %in% c("Intercept",
                                                                                            "intercept","precision","Precision"))]))
      if(resVar == "recruitment"){
        predicted=predicted/100
      }
    } 
    else {
      stop("Currently only ECCC 2011 and Johnson et al., 2020 models implemented")
    }
  }
  
  if (returnSample) {
    
    predLong = data.table::data.table(predictedTableSD)
    covTable$scnID = seq(1:nrow(covTable))
    predLong$scnID = seq(1:nrow(predLong))
    predLong = data.table::melt(predLong,
                    id.vars = "scnID",
                    variable.name = "replicate",
                    value.name = "value")
    
    resultDT = merge(covTable,predLong)    
    
  } 
  else {
    resultDT <- covTable
    resultDT$average = predicted
    resultDT$stdErr = predictedSD
    if (!is.null(nrow(qqs))) {
      resultDT$PIlow = qqs[,1]  
      resultDT$PIhigh = qqs[,2]  
    } else {
      resultDT$PIlow = qqs[[1]]  
      resultDT$PIhigh = qqs[[2]]
    }  
  }  
  tictoc::toc()
  
  return(resultDT)
}

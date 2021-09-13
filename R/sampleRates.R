#' Sample Demographic Rates
#'
#' Sample expected survival or recruitment rates based on samples of coefficient
#' values and optionally the model precision and interannual variation.
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
#' @param useQuantiles logical or numeric vector. Only relevant when
#'   \code{ignorePrecision = FALSE} and \code{returnSample = TRUE}. If
#'   \code{useQuantiles = TRUE}, each replicate population is assigned to a
#'   quantile of the distribution of variation around the expected values, and
#'   remains in that quantile as covariates change. If \code{useQuantiles} is a
#'   vector it is used as the quantiles. If \code{useQuantiles = FALSE},
#'   sampling is done independently for each combination of scenario and
#'   replicate, so the value for a particular replicate population in one
#'   scenario is unrelated to the values for that replicate in other scenarios.
#'   If interested in projecting impacts of changing disturbance on the
#'   trajectories of replicate populations set \code{useQuantiles = TRUE}.
#' @param interannualVar numeric or FALSE. Should the interannual variation in
#'   population growth be considered? If numeric it is the precision of the
#'   sample taken from the beta distribution.  TO DO: either remove this
#'   option from demographicParameters function, or update to align with
#'   popGrowthJohnson.
#' @param predInterval numeric vector with length 2. The prediction interval to
#'   use the default is 95% ie(\code{c(0.025,0.975)})
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
#'     \item{"PIlow"/"PIhigh"}{95% Prediction interval for estimated values}
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
                        useQuantiles,
                        interannualVar, 
                        predInterval = c(0.025,0.975)){
  
  tictoc::tic(paste0("Elapsed time for caribou prediction for ",
             resVar, " for ", modVer,":"))
  
  if ((!is.numeric(interannualVar) && interannualVar != FALSE) || 
      length(interannualVar) > 1) {
      stop("Expecting interannualVar to be a numeric precision parameter with length 1.")
    }
  
  
  whichCovariates <- names(coefValues)[!names(coefValues) %in% c("Intercept", 
                                                                   "intercept",
                                                                   "precision",
                                                                   "Precision")]
  
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
      if(sum(phi)==0){
        warning("Missing precision parameter.")
      }
      if(nrow(predictedTableSD)<1){
        stop("This code assumes at least one row.")
      }
      
      predictedTableSD = t(apply(predictedTableSD,1,betaSample,phi=phi,useQuantiles=useQuantiles))
    }
    
    if(interannualVar){
      pp = apply(predictedTableSD,1,betaSample,phi=interannualVar)
      
      pp=t(pp)
      predictedTableSD=pp
      
    }
    #Note: more relevant measure of uncertainty is perhaps bootstrap prediction interval
    #https://stats.stackexchange.com/questions/226565/bootstrap-prediction-interval
    #When precision is included points are not Gaussian distributed, and SD isn't an overly useful summary metric. 
    
    # Uncertainty across replicates
    predictedSD <- matrixStats::rowSds(predictedTableSD)
    
    if((length(useQuantiles)==1)&&!useQuantiles){
      qqs <- matrixStats::rowQuantiles(predictedTableSD,probs = predInterval)
    }else{
      if(length(setdiff(predInterval,useQuantiles))==0){
        colnames(predictedTableSD) <- useQuantiles
        
        selectCols = colnames(predictedTableSD)[grepl(predInterval[[1]],colnames(predictedTableSD),fixed=T)|grepl(predInterval[[2]],colnames(predictedTableSD),fixed=T)]
        if(length(selectCols)>length(predInterval)){
          stop("Handle case of replicates in useQuantiles.")
        }
        qqs <- subset(predictedTableSD,select=as.character(predInterval))
      }else{
        stop("If useQuantiles and predInterval are both specified, predInterval should be a subset of useQuantiles.")
      }
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
        if(sum(phi)==0){
          warning("Missing precision parameter.")
        }
        if(nrow(predictedTableSD)<1){
          stop("This code assumes at least one row.")
        }
        
        pp = apply(predictedTableSD,1,normalSample,sd=phi,useQuantiles=useQuantiles)
        
        predictedTableSD=t(pp)
        
      }
      
      # Uncertainty across replicates
      predictedSD <- matrixStats::rowSds(predictedTableSD)
      
      #TO DO: fix this - if using quantiles, pick appropriate quantiles from the set we have
      qqs <- matrixStats::rowQuantiles(predictedTableSD, probs = predInterval)
      
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
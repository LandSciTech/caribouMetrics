#' Create table of coefficients to use in population models, generated using bootstrapping
#'
#' @param covTable
#' @param coeffTable
#' @param coeffValues
#' @param modelType
#' @param model
#' @param ignorePrecision
#' @param returnSample
#' @param useQuantiles
#' @param interannualVar
#' 
#' @export

generatePopGrowthPredictions <- function(covTable,
                                         coeffTable,
                                         coeffValues,
                                         modelType,
                                         model,
                                         ignorePrecision,
                                         returnSample,
                                         useQuantiles,
                                         interannualVar){
  
  tic(paste0("Elapsed time for caribou prediction for ",
             model, " for ", modelType,":"))
  
  whichCovariates <- names(coeffValues)[!names(coeffValues) %in% c("Intercept", 
                                                                   "intercept",
                                                                   "precision",
                                                                   "Precision")]
  
  covTableRed <- covTable[, whichCovariates]
  
  if (grepl(x = modelType, pattern = "Johnson")) {
    intt <- coeffTable[, which(colnames(coeffTable) %in% c("Intercept", 
                                                           "intercept"))]
    phi <- coeffTable[, which(colnames(coeffTable) %in% c("Precision", 
                                                          "precision"))]
    
    predictedTableSD <- as.matrix(covTableRed) %*%
      t(coeffTable[,-which(colnames(coeffTable) %in% c("Intercept",
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
      
      pp = apply(predictedTableSD,1,betaSample,phi=phi,useQuantiles=useQuantiles)
      
      pp=t(pp)
      predictedTableSD=pp
      
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
    
    qqs <- matrixStats::rowQuantiles(predictedTableSD,probs = c(0.025,0.975))
    
    # Now the model calculations
    predicted <- as.numeric(exp(as.matrix(coeffValues)[
      which(colnames(coeffValues) %in% c("Intercept", "intercept"))] +
        as.matrix(covTableRed) %*% 
        as.matrix(coeffValues)[-which(colnames(coeffValues) %in% 
                                        c("Intercept", 
                                          "intercept",
                                          "precision",
                                          "Precision"))]))
  } 
  else {
    if (grepl(x = modelType, pattern = "ECCC")) {
      
      intt <- coeffTable[, which(colnames(coeffTable) %in% c("Intercept", "intercept"))]
      phi <- coeffTable[, which(colnames(coeffTable) %in% c("Precision", 
                                                            "precision"))]
      
      predictedTableSD <- as.matrix(covTableRed) %*%
        t(coeffTable[,-which(colnames(coeffTable) %in% c("Intercept","intercept","precision","Precision"))])
      
      predictedTableSD<-t(apply(predictedTableSD,1,function(x,intt){(intt+x)},intt=intt))
      
      if(model=="recruitment"){
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
        
        pp=t(pp)
        predictedTableSD=pp
        
      }
      
      # Uncertainty across replicates
      predictedSD <- matrixStats::rowSds(predictedTableSD)
      
      qqs <- matrixStats::rowQuantiles(predictedTableSD,probs=c(0.025,0.975))
      
      # Now the model calculations
      predicted <- as.numeric((as.matrix(coeffValues)[which(colnames(coeffValues) %in% c("Intercept",
                                                                                         "intercept"))] +
                                 as.matrix(covTableRed) %*%
                                 as.matrix(coeffValues)[-which(colnames(coeffValues) %in% c("Intercept",
                                                                                            "intercept","precision","Precision"))]))
      if(model=="recruitment"){
        predicted=predicted/100
      }
    } 
    else {
      stop("Currently only ECCC 2011 and Johnson et al., 2020 models implemented")
    }
  }
  
  if (returnSample) {
    
    predLong = data.table(predictedTableSD)
    covTable$scnID = seq(1:nrow(covTable))
    predLong$scnID = seq(1:nrow(predLong))
    predLong = melt(predLong,
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
    } 
    else {
      resultDT$PIlow = qqs[[1]]  
      resultDT$PIhigh = qqs[[2]]
    }  
  }  
  toc()
  
  return(resultDT)
}

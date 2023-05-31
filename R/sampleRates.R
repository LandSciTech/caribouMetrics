# Copyright 2021 Tati Micheletti & Her Majesty the Queen in Right of Canada as represented by the Minister of the Environment
# License GPL-3
#NOTICE: This function has been modified from https://github.com/tati-micheletti/caribouPopGrowthModel/blob/master/R/generatePopGrowthPredictions.R

#' Sample demographic rates
#'
#' 
#' @param coefSamples matrix. Bootstrapped coefficients with one row per
#'   replicate and one column per coefficient
#' @param coefValues data.table. One row table with expected values for each
#'   coefficient
#' @inheritParams demographicCoefficients
#' @param quantilesToUse numeric vector of length `coefSamples`. See `useQuantiles`.
#' @param transformFn function used to transform demographic rates.
#'
#' @inheritParams demographicRates
#'
#' @return For `sampleRates` a similar data frame for one response variable
#'
#' @examples 
#' cfs <- getCoefs(popGrowthTableJohnsonECCC, "recruitment", "Johnson", "M3")
#' 
#' cfSamps <- sampleCoefs(cfs[[1]], 10)
#'
#' # disturbance scenarios
#' distScen <- data.frame(Total_dist = 1:10/10)
#'
#' # return summary across replicates
#' sampleRates(distScen, cfSamps$coefSamples, cfSamps$coefValues,
#'             "Johnson", "recruitment", ignorePrecision = TRUE, 
#'             returnSample = FALSE)
#'
#' # return one row per replicate * scenario
#' sampleRates(distScen, cfSamps$coefSamples, cfSamps$coefValues,
#'             "Johnson", "recruitment", ignorePrecision = TRUE, 
#'             returnSample = TRUE)
#'
#' # return one row per replicate * scenario with replicates assigned to a quantile
#' sampleRates(distScen, cfSamps$coefSamples, cfSamps$coefValues,
#'             "Johnson", "recruitment", ignorePrecision = TRUE, 
#'             returnSample = TRUE, 
#'             quantilesToUse = quantile(x = c(0, 1),
#'                                       probs = seq(0.025, 0.975, length.out = 10)))
#'
#'
#' @rdname demographicRates
#' @export

sampleRates <- function(covTable,
                        coefSamples,
                        coefValues,
                        modelVersion,
                        resVar,
                        ignorePrecision,
                        returnSample,
                        quantilesToUse=NULL,
                        predInterval = c(0.025,0.975),
                        transformFn = function(y){y}){

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

  # set up components that are different for different model versions
  if (grepl(x = modelVersion, pattern = "Johnson")) {
    predictSDFun <- function(x, intt){transformFn(exp(intt + x))}
    phiSampleFun <- betaSample
    predictFun <- function(coefValues, covTableRed){transformFn(
      exp(as.numeric(as.matrix(coefValues)[
        which(colnames(coefValues) %in% c("Intercept", "intercept"))] +
          as.matrix(covTableRed) %*%
          as.matrix(coefValues)[-which(colnames(coefValues) %in%
                                         c("Intercept",
                                           "intercept",
                                           "precision",
                                           "Precision"))]))
    )}
    recruitDiv <- 1
  } else if (grepl(x = modelVersion, pattern = "ECCC")) {
    predictSDFun <- function(x, intt){intt + x}
    phiSampleFun <- normalSample
    predictFun <- function(coefValues, covTableRed){
      as.numeric(as.matrix(coefValues)[
        which(colnames(coefValues) %in% c("Intercept", "intercept"))] +
          as.matrix(covTableRed) %*%
          as.matrix(coefValues)[-which(colnames(coefValues) %in%
                                         c("Intercept",
                                           "intercept",
                                           "precision",
                                           "Precision"))])
    }
    recruitDiv <- 100
  }else {
    stop("Currently only ECCC 2011 and Johnson et al., 2020 models implemented")
  }

  intt <- coefSamples[, which(colnames(coefSamples) %in% c("Intercept",
                                                           "intercept"))]
  phi <- coefSamples[, which(colnames(coefSamples) %in% c("Precision",
                                                          "precision"))]

  predictedTableSD <- as.matrix(covTableRed) %*%
    t(coefSamples[,-which(colnames(coefSamples) %in% c("Intercept",
                                                       "intercept",
                                                       "precision",
                                                       "Precision"))])

  predictedTableSD <- t(apply(predictedTableSD, 1, predictSDFun, intt = intt))

  if (resVar == "recruitment") {
    predictedTableSD <- predictedTableSD/recruitDiv
  }
  #predictedTableSD is bootstrap sample of mean.
  #Precision parameter gives info about the amount of scatter around the mean.
  #Easier to understand the principles by thinking about a simpler linear regression example.
  #Using the notation in the regression analysis section of this wikipedia article, predictedTableSD is sample from distribution of yhat_d, rather than y_d: https://en.wikipedia.org/wiki/Prediction_interval
  if(!ignorePrecision){
    if(length(phi) == 0){
      stop("Missing precision parameter. Set ignorePrecision = TRUE or",
           " add a precision column to coefSamples", call. = FALSE)
    }
    if(nrow(predictedTableSD)<1){
      stop("This code assumes at least one row.")
    }

    predictedTableSD <- t(apply(predictedTableSD, 1, phiSampleFun, phi = phi,
                                quantilesToUse = quantilesToUse))
  }

  #Note: more relevant measure of uncertainty is perhaps bootstrap prediction interval
  #https://stats.stackexchange.com/questions/226565/bootstrap-prediction-interval
  #When precision is included points are not Gaussian distributed, and SD isn't an overly useful summary metric.

  # Uncertainty across replicates
  predictedSD <- matrixStats::rowSds(predictedTableSD)

  if(is.null(quantilesToUse)){
    qqs <- matrixStats::rowQuantiles(predictedTableSD, probs = predInterval)
  }else{
    qqs <- matrixStats::rowQuantiles(predictedTableSD, probs = c(0,1))
  }
  # Now the model calculations
  predicted <- predictFun(coefValues, covTableRed)

  if(resVar == "recruitment"){
    predicted = predicted/recruitDiv
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

  return(resultDT)
}

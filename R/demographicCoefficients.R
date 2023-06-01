#' Sample demographic regression model coefficients
#'
#'
#' Select the regression coefficient values and standard errors for the desired
#' model version (see `popGrowthTableJohnsonECCC` for options) and then sample
#' from the Gaussian distribution for each replicate population.
#' `demographicCoefficients` is a wrapper around [getCoefs()], which selects
#' coefficients and [sampleCoefs()], which samples coefficients, for both the
#' survival and recruitment models.
#'
#' Each population is optionally assigned to quantiles of the error
#' distributions for survival and recruitment. Using quantiles means that the
#' population will stay in these quantiles as disturbance changes over time, so
#' there is persistent variation in recruitment and survival among example
#' populations. See [demographicRates()] for more details. 
#'
#' @param replicates integer. Number of replicate populations.
#' @param modelVersion character. Which model version to use. Currently the only
#'   option is "Johnson" for the model used in Johnson et. al. (2020), but
#'   additional options may be added in the future.
#' @param survivalModelNumber,recruitmentModelNumber character. Which model
#'   number to use see [popGrowthTableJohnsonECCC] for options.
#' @param useQuantiles logical or numeric. If it is a numeric vector it must be
#'   length 2 and give the low and high limits of the quantiles to use. If
#'   `useQuantiles != FALSE`, each replicate population is assigned to a
#'   quantile of the distribution of variation around the expected values, and
#'   remains in that quantile as covariates change. If `useQuantiles = TRUE`,
#'   replicate populations will be assigned to quantiles in the default range of
#'   0.025 and 0.975.  
#' @param populationGrowthTable data.frame.[popGrowthTableJohnsonECCC] is
#'   included in the package and should be used in most cases. A custom table of
#'   model coefficients and standard errors or confidence intervals can be
#'   provided but it must match the column names of [popGrowthTableJohnsonECCC].
#'   If the table does not contain the standard error it is calculated from the
#'   confidence interval.
#'
#' @return For `demographicCoefficients` a list with elements:
#'   * "modelVersion": The name of the model version
#'   * "coefSamples_Survival" and"coefSamples_Recruitment":
#'   lists with elements:
#'     * "coefSamples": Bootstrapped coefficients with `replicates` rows
#'     * "coefValues": Coefficient values taken from `populationGrowthTable`
#'     * "quantiles": A vector of randomly selected quantiles between 0.025 and
#'   0.975 with length `replicates`
#'
#' @references Johnson, C.A., Sutherland, G.D., Neave, E., Leblond, M., Kirby,
#'   P., Superbie, C. and McLoughlin, P.D., 2020. Science to inform policy:
#'   linking population dynamics to habitat for a threatened species in Canada.
#'   Journal of Applied Ecology, 57(7), pp.1314-1327.
#'   <https://doi.org/10.1111/1365-2664.13637>
#'
#' @examples
#' # sample coefficients for default models
#' demographicCoefficients(10)
#'
#' # try a different model
#' demographicCoefficients(10, modelVersion = "Johnson", survivalModelNumber = "M1",
#'                         recruitmentModelNumber = "M3")
#'
#' @family demography
#' @export

demographicCoefficients <- function(replicates,
                           modelVersion = "Johnson",
                           survivalModelNumber = "M1",
                           recruitmentModelNumber = "M4",
                           useQuantiles = TRUE,
                           populationGrowthTable = popGrowthTableJohnsonECCC){

  if(!all(colnames(popGrowthTableJohnsonECCC) %in%
          colnames(populationGrowthTable))) {
    stop("populationGrowthTable must contain all colnames in popGrowthTableJohnsonECCC")
  }

  if(any(length(modelVersion) > 1, length(survivalModelNumber) > 1,
         length(recruitmentModelNumber) > 1)){
    stop("Multiple models. modelVersion, survivalModelNumber, ",
         "and recruitmentModelNumber must have length 1", call. = FALSE)
  }

  quantsToUse <- prepQuantiles(useQuantiles)

  populationGrowthTable <- data.table::data.table(populationGrowthTable)

  DT_S <- getCoefs(populationGrowthTable,
                    resVar = "femaleSurvival",
                    modelVersion = modelVersion,
                    modNum = survivalModelNumber)[[1]]

  coefSamples_S <- sampleCoefs(DT_S, replicates)

  DT_R <- getCoefs(populationGrowthTable,
                    resVar = "recruitment",
                    modelVersion = modelVersion,
                    modNum = recruitmentModelNumber)[[1]]

  coefSamples_R <- sampleCoefs(DT_R, replicates)

  if(!is.null(quantsToUse)){
  coefSamples_S$quantiles <- sample(getQuantiles(nrow(coefSamples_S$coefSamples),
                                                 low = quantsToUse[1],
                                                 high = quantsToUse[2]),
                                    replace = FALSE)
  coefSamples_R$quantiles <- sample(getQuantiles(nrow(coefSamples_R$coefSamples),
                                                 low = quantsToUse[1],
                                                 high = quantsToUse[2]),
                                    replace = FALSE)
}

  return(list(modelVersion = modelVersion,
              coefSamples_Survival = coefSamples_S,
              coefSamples_Recruitment = coefSamples_R))
}

getQuantiles<-function(x,low=0.025,high=0.975){
  return(low+(seq(0,x-1)/(x-1))*(high-low))
}

prepQuantiles <- function(useQuantiles, quantilesIn = NULL){
  if(length(useQuantiles) == 2){
    if(!all(min(useQuantiles) >= 0, max(useQuantiles) <= 1)){
      stop("useQuantiles must be between 0 and 1")
    }
    if(!is.null(quantilesIn)){
      warning("popGrowthPars contains quantiles so they are used and useQuantiles is ignored",
              call. = FALSE)
      
      quantsToUse <- quantilesIn
    }else{
      quantsToUse <- useQuantiles
    }
  } else if(length(useQuantiles) == 1){
    if(useQuantiles){
      if(!is.null(quantilesIn)){
        message("popGrowthPars contains quantiles so they are used instead of the defaults")
        
        quantsToUse <- quantilesIn
      }else{
        quantsToUse <- c(0.025, 0.975)
      }
    } else {
      quantsToUse <- NULL
    }
  } else {
    stop("useQuantiles must have length 1 or 2")
  }

  return(quantsToUse)
}

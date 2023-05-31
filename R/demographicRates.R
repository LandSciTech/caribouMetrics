#' Sample expected survival and recruitment rates
#'
#' @description 
#' Apply sampled coefficients to disturbance covariates to calculate expected
#' recruitment (\eqn{\bar{R}_t}) and survival (\eqn{\bar{S}_t}) according to the
#' beta regression models estimated by Johnson et al. (2020): \deqn{\bar{R}_t
#' \sim Beta(\mu^R_t,\phi^R);
#' log(\mu^R_t)=\dot{\beta^R_0}+\dot{\beta^R_a}A_t+\dot{\beta}^R_fF_t}
#'
#' \deqn{\bar{S}_t \sim (46\times
#' Beta(\mu^S_t,\phi^S)-0.5)/45;log(\mu^S_t)=\dot{\beta^S_0}+\dot{\beta^S_a}A_t}
#'
#' \eqn{\phi^R} and \eqn{\phi^S} are the precisions of the Beta distributed
#' errors and are only available for some model versions, or can optionally be
#' ignored.
#'
#' `demographicRates` is a wrapper around `sampleRates` to sample both survival
#' and recruitment rates based on the result of [demographicCoefficients()] and
#' using recommended defaults.
#' 
#' @details 
#' 
#' TODO: confirm this is the right equation
#' ## Recruitment
#' Recruitment is modeled as
#' \deqn{\log(R_t)=\beta^R_0+\beta^R_a A_t+\beta^R_f F_t+\epsilon^R_t; \epsilon^R_t \sim \text{Normal}(0,\sigma^2_{R})}
#' 
#' ## Survival
#' \deqn{S_t=\text{min}(46 \tilde{S_t}-0.5)/45,1); \log(\tilde{S_t})=\beta^S_0+\beta^S_a A_t+\epsilon^S_t; \epsilon^S_t \sim \text{Normal}(0,\sigma^2_{S})}
#'
#' 
#' @param covTable data.frame. A table of covariate values to be used. Column
#'   names must match the coefficient names in [popGrowthTableJohnsonECCC]. Each
#'   row is a different scenario.
#' @param popGrowthPars list. Coefficient values and (optionally) quantiles
#'   returned by `demographicCoefficients`.
#' @param ignorePrecision logical. Should the precision of the model be used if
#'   it is available? When precision is used variation among populations around
#'   the National mean responses is considered in addition to the uncertainty
#'   about the coefficient estimates.
#' @param returnSample logical. If TRUE the returned data.frame has replicates *
#'   scenarios rows. If FALSE the returned data.frame has one row per scenario
#'   and additional columns summarizing the variation among replicates. See
#'   Value for details.
#' @param useQuantiles logical or numeric. If it is a numeric vector it must be
#'   length 2 and give the low and high limits of the quantiles to use. Only
#'   relevant when `ignorePrecision = FALSE`. If `useQuantiles != FALSE`, each
#'   replicate population is assigned to a quantile of the distribution of
#'   variation around the expected values, and remains in that quantile as
#'   covariates change. If `useQuantiles != FALSE` and popGrowthPars contains
#'   quantiles, those quantiles will be used. If `useQuantiles = TRUE` and
#'   popGrowthPars does not contain quantiles, replicate populations will be
#'   assigned to quantiles in the default range of 0.025 and 0.975. If
#'   `useQuantiles = FALSE`, sampling is done independently for each combination
#'   of scenario and replicate, so the value for a particular replicate
#'   population in one scenario is unrelated to the values for that replicate in
#'   other scenarios. Useful for projecting impacts of changing disturbance on
#'   the trajectories of replicate populations.
#' @param predInterval numeric vector with length 2. The default 95% interval
#'   is (`c(0.025,0.975)`). Only relevant when `returnSample = TRUE` and
#'   `quantilesToUse = NULL`.
#' @param transformFns list of functions used to transform demographic rates.
#'   The default is `list(S_transform = function(y){(y*46-0.5)/45},R_transform =
#'   function(y){y})`. The back transformation is applied to survival rates as
#'   in Johnson et al. 2020.
#'
#' @return A data.frame of predictions. The data.frame includes all columns in
#'   `covTable` with additional columns depending on `returnSample`.
#'
#'   If `returnSample = FALSE` the number of rows is the same as the number of
#'   rows in `covTable`, additional columns are:
#'   * "S_bar" and "R_bar": The mean estimated values of survival and
#'   recruitment (calves per cow)
#'   * "S_stdErr" and "R_stdErr": Standard error of the estimated values
#'   * "S_PIlow"/"S_PIhigh" and "R_PIlow"/"R_PIhigh": If not using quantiles,
#'   95\% of values fall within this range. If using quantiles, maximum and
#'   minimum values are returned.
#'
#'   If `returnSample = TRUE` the number of rows is `nrow(covTable) *
#'   replicates` additional columns are:
#'   * "scnID": A unique identifier for scenarios provided in
#'   `covTable`
#'   * "replicate": A replicate identifier, unique within each scenario
#'   * "S_bar" and "R_bar": The expected values of survival and
#'   recruitment (calves per cow)
#'
#'
#' @examples
#' # get coefficient samples
#' coefs <- demographicCoefficients(10)
#'
#' # table of different scenarios to test
#' covTableSim <- expand.grid(Anthro = seq(0, 90, by = 20),
#'                            fire_excl_anthro = seq(0, 70, by = 20))
#' covTableSim$Total_dist = covTableSim$Anthro + covTableSim$fire_excl_anthro
#'
#' demographicRates(covTableSim, coefs)
#'
#' @family demography
#' @export
demographicRates <- function(covTable,
                             popGrowthPars,
                             ignorePrecision = FALSE,
                             returnSample = FALSE,
                             useQuantiles = TRUE,
                             predInterval = list(PI_R = c(0.025,0.975),
                                                 PI_S = c(0.025,0.975)),
                             transformFns = list(S_transform = function(y){(y*46-0.5)/45},R_transform = function(y){y})){

  quantsToUse <- prepQuantiles(useQuantiles, popGrowthPars$coefSamples_Survival$quantiles)
  if(is.null(quantsToUse)){
    popGrowthPars$coefSamples_Survival$quantiles <- NULL
  }else{
    q <- getQuantiles(nrow(popGrowthPars$coefSamples_Survival$coefSamples),
                      low = quantsToUse[1],
                      high = quantsToUse[2])

    if(is.null(popGrowthPars$coefSamples_Survival$quantiles)){
      popGrowthPars$coefSamples_Survival$quantiles <- sample(q, replace = FALSE)
    }
  }

  quantsToUse <- prepQuantiles(useQuantiles, popGrowthPars$coefSamples_Recruitment$quantiles)
  if(is.null(quantsToUse)){
    popGrowthPars$coefSamples_Recruitment$quantiles <- NULL
  }else{
    q <- getQuantiles(nrow(popGrowthPars$coefSamples_Recruitment$coefSamples),
                      low = quantsToUse[1],
                      high = quantsToUse[2])
    if(is.null(popGrowthPars$coefSamples_Recruitment$quantiles)){
      popGrowthPars$coefSamples_Recruitment$quantiles <- sample(q, replace = FALSE)
    }
  }

  pred_S <- sampleRates(covTable = covTable,
                        coefSamples = popGrowthPars$coefSamples_Survival[["coefSamples"]],
                        coefValues = popGrowthPars$coefSamples_Survival[["coefValues"]],
                        modelVersion =  popGrowthPars$modelVersion,
                        quantilesToUse = popGrowthPars$coefSamples_Survival[["quantiles"]],
                        resVar = "femaleSurvival",
                        ignorePrecision = ignorePrecision,
                        returnSample = returnSample,
                        predInterval = predInterval[["PI_S"]],
                        transformFn = transformFns$S_transform)
  pred_R <- sampleRates(covTable = covTable,
                        coefSamples = popGrowthPars$coefSamples_Recruitment[["coefSamples"]],
                        coefValues = popGrowthPars$coefSamples_Recruitment[["coefValues"]],
                        quantilesToUse = popGrowthPars$coefSamples_Recruitment[["quantiles"]],
                        modelVersion =  popGrowthPars$modelVersion,
                        resVar = "recruitment",
                        ignorePrecision = ignorePrecision,
                        returnSample = returnSample,
                        predInterval = predInterval[["PI_R"]],
                        transformFn = transformFns$R_transform)
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

  if(returnSample){
    rateSamples<-rateSamples[order(rateSamples$replicate, rateSamples$scnID),]
  }

  return(rateSamples)
}

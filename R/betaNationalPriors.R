#' Get prior parameters for Bayesian Beta demographic rate models 
#'
#' Returns prior parameter values for Bayesian Beta demographic rate models described by [Hughes et al. (2025)](https://doi.org/10.1016/j.ecoinf.2025.103095). 
#' The starting point is estimated coefficients from national
#' demographic-disturbance relationships in the table
#' `popGrowthTableJohnsonECCC`. 
#' 
#' Standard errors and random effects of year have
#' been calibrated so that the 95% prior prediction intervals for survival and
#' recruitment from the Bayesian model match the range between the 2.5% and
#' 97.5% quantiles of 1000 survival and recruitment trajectories from the
#' national demographic model [caribouPopGrowth()]. A prior for the
#' unknown composition survey bias correction term `c` is set by specifying an
#' apparent number of adult females per collared animal(`cowMult`) and minimum
#' and maximum values for each of the ratio of bulls to cows (\eqn{q}), the
#' probability of misidentifying young bulls as adult females and vice versa
#' (\eqn{u}), and the probability of missing calves (\eqn{z}) in composition
#' surveys. See [compositionBiasCorrection()] and [Hughes et al. (2025)](https://doi.org/10.1016/j.ecoinf.2025.103095) for additional details.
#' 
#' @param modList a named list of modifiers to use to change the priors. If a
#'   modifier is supplied here the corresponding argument below is ignored.
#' @param returnValues logical. Default is TRUE. If FALSE returns strings for
#'   some values showing the initial values and the modifier ie "0.9 * 1.05"
#' @param sAnthroSlopeSE Standard deviation of effect of disturbance on survival.
#' @param rAnthroSlopeSE Standard deviation of effect of anthropogenic disturbance on recruitment.
#' @param rFireSlopeSE Standard deviation of effect of fire on recruitment.
#' @param sIntSE Standard deviation of survival intercept.
#' @param rIntSE Standard deviation of recruitment intercept.
#' @param sNuMin,sNuMax Uniform prior for coefficient of variation among years for recruitment. 
#' @param rNuMin,rNuMax Uniform prior for coefficient of variation among years for recruitment. 
#' @param qMin number in 0, 1. Minimum ratio of bulls to cows in composition
#'   survey groups.
#' @param qMax number in 0, 1. Maximum ratio of bulls to cows in composition
#'   survey groups.
#' @param uMin number in 0, 1. Minimum probability of misidentifying young bulls
#'   as adult females and vice versa in composition survey.
#' @param uMax number in 0, 1. Maximum probability of misidentifying young bulls
#'   as adult females and vice versa in composition survey.
#' @param zMin number in 0, 1. Minimum probability of missing calves in
#'   composition survey.
#' @param zMax number in 0, 1. Maximum probability of missing calves in
#'   composition survey.
#' @param cowMult number. The apparent number of adult females per collared
#'   animal in composition survey.
#' @inheritParams subsetNationalCoefs
#' @inheritParams getNationalCoefficients
#'
#' @return a list with values:
#'
#' * R_b0_mu: Recruitment intercept
#' * R_b0_sd: Recruitment intercept standard deviation,
#' * R_b1_mu: Recruitment anthropogenic disturbance slope,
#' * R_b1_sd: Recruitment anthropogenic disturbance standard deviation,
#' * R_b2_mu: Recruitment fire excluding anthropogenic disturbance
#'   slope,
#' * R_b2_sd: Recruitment fire excluding anthropogenic disturbance standard deviation,
#' * R_cv_min: Min of the prior distribution of the random effect of year
#'   on recruitment,
#' * R_cv_max: Max of the prior distribution of the random effect of year on recruitment,
#' * S_b0_mu: Adult female survival intercept,
#' * S_b0_sd: Adult female survival intercept standard error times modifier,
#' * S_b1_mu: Adult female survival anthropogenic disturbance slope,
#' * S_b1_sd: Adult female survival anthropogenic disturbance standard deviation,
#' * S_cv_min: Min of the prior distribution of the random effect of year
#'   on adult female survival,
#' * S_cv_max: Max of the prior distribution of the random
#'   effect of year on adult female survival,
#' * qMin,qMax,uMin,uMax,zMin,zMax,cowMult: Composition bias parameters.   
#' 
#' @references    
#'   Hughes, J., Endicott, S., Calvert, A.M. and Johnson, C.A., 2025.
#'   Integration of national demographic-disturbance relationships and local
#'   data can improve caribou population viability projections and inform
#'   monitoring decisions. Ecological Informatics, 87, p.103095.
#'   <https://doi.org/10.1016/j.ecoinf.2025.103095>
#'   
#' @examples
#' betaNationalPriors()
#' 
#' @seealso [bbouNationalPriors::bbouNationalPriors()] for priors for bboutools
#'   models from national demographic-disturbance relationships.
#'
#' @family demography
#' @export
betaNationalPriors <- function(modList = NULL,
                      survivalModelNumber = "M1",
                      recruitmentModelNumber = "M4",
                      rAnthroSlopeSE = 0.006,
                      rFireSlopeSE = 0.002,
                      sAnthroSlopeSE = 0.0005,
                      sIntSE = 0.06,
                      sNuMin = 0.01,
                      sNuMax = 0.13,
                      rIntSE = 0.35,
                      rNuMin =0.01,
                      rNuMax =0.7,
                      qMin=0, qMax =0, 
                      uMin = 0, uMax = 0, 
                      zMin = 0, zMax = 0, 
                      cowMult = 6,
                      populationGrowthTable = caribouMetrics::popGrowthTableJohnsonECCC,
                      modelVersion = "Johnson",
                      r.inv.link="exp",
                      returnValues = TRUE) {
  # modList=paramTable

  expectMods <- c(as.list(environment()))
  expectMods$modList <- NULL
  if (is.null(modList)) {
    modList <- expectMods
  } else {
    # keep all values in modifiers and add any that are missing using values in
    # expectMods
    modList <- c(modList, expectMods[which(!names(expectMods) %in% names(modList))])
  }

  popGrowthPars <- getNationalCoefficients(
    2,
    modelVersion = modelVersion,
    survivalModelNumber = modList$survivalModelNumber,
    recruitmentModelNumber = modList$recruitmentModelNumber,
    populationGrowthTable = populationGrowthTable
  )

  rPriorCoefs <- popGrowthPars$coefSamples_Recruitment$coefValues
  rPriorStdErrs <- popGrowthPars$coefSamples_Recruitment$coefStdErrs
  sPriorCoefs <- popGrowthPars$coefSamples_Survival$coefValues
  sPriorStdErrs <- popGrowthPars$coefSamples_Survival$coefStdErrs

  # check for variables not programmed in jags
  diff_pars <- setdiff(
    names(rPriorCoefs),
    c("Intercept", "Precision", "Anthro", "fire_excl_anthro")
  )
  if (length(diff_pars) > 0) {
    stop(
      "The recruitment model contains unrecognized coefficients: ",
      paste0(diff_pars, collapse = ", ")
    )
  }

  diff_pars <- setdiff(
    names(sPriorCoefs),
    c("Intercept", "Precision", "Anthro", "fire_excl_anthro")
  )
  if (length(diff_pars) > 0) {
    stop(
      "The survival model contains unrecognized coefficients: ",
      paste0(diff_pars, collapse = ", ")
    )
  }
  
  #####
  #get bias coefficient priors - lognormally distributed
  #nr=10000
  #cs = compositionBiasCorrection(w=modList$cowMult,q=runif(nr,modList$qMin,modList$qMax),
  #                              u=runif(nr,modList$uMin,modList$uMax),
  #                              z=runif(nr,modList$zMin,modList$zMax),approx=T)
  #bias.Prior1 = cs$mu
  #bias.Prior2 = cs$sig2^0.5
  compositionBiasPars <- modList[c("cowMult","qMin","qMax","uMin","uMax","zMin","zMax")]

  if (returnValues) {
    betaPriors <- list(
      R_b0_mu = rPriorCoefs$Intercept,
      R_b0_sd = modList$rIntSE,
      R_b1_mu = rPriorCoefs$Anthro,
      R_b1_sd = modList$rAnthroSlopeSE,
      R_b2_mu = rPriorCoefs$fire_excl_anthro,
      R_b2_sd = modList$rFireSlopeSE,
      R_cv_min = modList$rNuMin,
      R_cv_max = modList$rNuMax,
      S_b0_mu = sPriorCoefs$Intercept,
      S_b0_sd = modList$sIntSE,
      S_b1_mu = sPriorCoefs$Anthro,
      S_b1_sd = modList$sAnthroSlopeSE,
      S_cv_min = modList$sNuMin,
      S_cv_max = modList$sNuMax,
      R_inv_link =modList$r.inv.link
    )

    # replace NULL values with 0 for when anthro or fire is not included
    betaPriors <- lapply(betaPriors, function(x) {
      if (is.null(x) || length(x) == 0) {
        1e-10
      } else {
        x
      }
    })
    betaPriors <- c(betaPriors,compositionBiasPars)
  } else {
    betaPriors <- list(
      R_b0_mu = rPriorCoefs$Intercept,
      R_b0_sd = round(modList$rIntSE, 4),
      R_b1_mu = rPriorCoefs$Anthro,
      R_b1_sd = round(modList$rAnthroSlopeSE, 4),
      R_b2_mu = rPriorCoefs$fire_excl_anthro,
      R_b2_sd = round(modList$rFireSlopeSE,4),
      R_cv_min = modList$rNuMin,
      R_cv_max = modList$rNuMax,
      S_b0_mu = sPriorCoefs$Intercept,
      S_b0_sd = round(modList$sIntSE, 4),
      S_b1_mu = sPriorCoefs$Anthro,
      S_b1_sd = round(modList$sAnthroSlopeSE, 4),
      S_cv_min = modList$sNuMin,
      S_cv_max = modList$sNuMax
    )
  }
  
  return(betaPriors)
}

simCovariates <- function(initAnthro, initFire, numYears, anthroSlope,
                          anthroSlopeFuture, futureStep, fireSlope = 0) {
  
  iv <- c(initAnthro,initFire,anthroSlope,anthroSlopeFuture,fireSlope)
  iv <- iv[!is.na(iv)]
  if(length(iv)!=5){
    covariates <- data.frame(time=NA)
    covariates$Anthro <- NA
    covariates$fire_excl_anthro <- NA
    covariates$Total_dist <- NA
    covariates <- subset(covariates,!is.na(time))
    return(covariates)
  }
  covInit <- data.frame(Anthro = initAnthro, fire_excl_anthro = initFire)
  for (t in 1:numYears) {
    # t=1s
    if (t == 1) {
      cov <- covInit
    } else {
      cov <- covPrev
    }
    if (t > 1) {
      if (t >= futureStep) {
        cov$Anthro <- covPrev$Anthro + anthroSlopeFuture
      } else {
        cov$Anthro <- covPrev$Anthro + anthroSlope
      }
    }
    covPrev <- cov
    if (fireSlope == 0) {
      cov$fire_excl_anthro <- pmax(0, rnorm(1, cov$fire_excl_anthro, 0.0001))
    } else {
      cov$fire_excl_anthro <- cov$fire_excl_anthro + fireSlope * (t - 1)
    }
    cov$time <- t
    if (t == 1) {
      covariates <- cov
    } else {
      covariates <- rbind(covariates, cov)
    }
  }
  covariates$Anthro <- pmin(100, pmax(0, covariates$Anthro))
  covariates$fire_excl_anthro <- pmin(100, pmax(0, covariates$fire_excl_anthro))
  covariates$Total_dist <- pmin(100, covariates$Anthro + covariates$fire_excl_anthro)

  return(covariates)
}


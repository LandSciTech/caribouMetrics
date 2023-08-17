#' Get model priors from national models
#'
#'
#'
#' @param modList a named list of modifiers to use to change the priors. If a
#'   modifier is supplied here the corresponding argument below is ignored.
#' @param returnValues logical. Default is TRUE. If FALSE returns strings for
#'   some values showing the initial values and the modifier ie "0.9 * 1.05"
#' @param sAnthroSlopeSEMod anthropogenic disturbance slope survival uncertainty multiplier. 
#'   1 - 10
#' @param rAnthroSlopeSEMod anthropogenic disturbance slope recruitment uncertainty
#'   multiplier. 1 - 10
#' @param sIntSEMod survival intercept uncertainty multiplier. 1 - 10
#' @param rIntSEMod recruitment intercept uncertainty multiplier. 1 - 10
#' @param sInterannualVar interannual coefficient of variation for survival. 0-1. See
#'   [caribouPopGrowth()] and functions therein for details
#' @param sInterannualVarSE uncertainty about interannual variation in survival. 0-1
#' @param rInterannualVar interannual coefficient of variation for recruitment. 0-1
#' @param rInterannualVarSE uncertainty about interannual variation in recruitment. 0-1.
#' @param qMin number in 0, 1. Minimum ratio of bulls to cows in composition survey groups.
#' @param qMax number in 0, 1. Maximum ratio of bulls to cows in composition survey groups.
#' @param uMin number in 0, 1. Minimum probability of misidentifying young bulls as adult females and vice versa in composition survey.
#' @param uMax number in 0, 1. Maximum probability of misidentifying young bulls as adult females and vice versa in composition survey.
#' @param zMin number in 0, 1. Minimum probability of missing calves in composition survey.
#' @param zMax number in 0, 1. Maximum probability of missing calves in composition survey.
#' @param cowMult number. The apparent number of adult females per collared animal in composition survey.
#' @inheritParams getCoefs
#' @inheritParams demographicCoefficients
#'
#' @return a list with values:
#' 
#' * l.R.Prior1: Recruitment intercept  
#' * l.R.Prior2: Recruitment intercept standard error times modifier,
#' * beta.Rec.anthro.Prior1: Recruitment anthropogenic disturbance slope,
#' * beta.Rec.anthro.Prior2: Recruitment anthropogenic disturbance standard
#'   error times modifier,
#' * beta.Rec.fire.Prior1: Recruitment fire excluding anthropogenic disturbance
#'   slope,
#' * beta.Rec.fire.Prior2: Recruitment fire excluding anthropogenic disturbance
#'   standard error,
#' * sig.R.Prior1: Interannual coefficient of variation for recruitment,
#' * sig.R.Prior2: uncertainty about interannual variation in recruitment,
#' * l.Saf.Prior1: Adult female survival intercept,
#' * l.Saf.Prior2: Adult female survival intercept standard error times modifier,
#' * beta.Saf.Prior1: Adult female survival anthropogenic disturbance slope,
#' * beta.Saf.Prior2: Adult female survival anthropogenic disturbance standard
#'   error times modifier,
#' * sig.Saf.Prior1: Interannual coefficient of variation for adult female survival,
#' * sig.Saf.Prior2: Uncertainty about interannual variation in adult female survival,
#' * bias.Prior1: Mean composition survey bias correction term,
#' * bias.Prior2: Standard deviation of composition survey bias correction term
#' 
#' @examples
#' getPriors()
#' 
#' @family demography
#' @export
getPriors <- function(modList = NULL,
                      survivalModelNumber = "M1",
                      recruitmentModelNumber = "M4",
                      rAnthroSlopeSEMod = 4,
                      sAnthroSlopeSEMod = 5,
                      sIntSEMod = 5,
                      sInterannualVar = 0.08696 * 0.4,
                      sInterannualVarSE = 0.03,
                      rIntSEMod = 3,
                      rInterannualVar = 0.46 * 0.5,
                      rInterannualVarSE = 0.22,
                      qMin=0, qMax =0.6, 
                      uMin = 0, uMax = 0.2, 
                      zMin = 0, zMax = 0.2, 
                      cowMult = 6,
                      populationGrowthTable = caribouMetrics::popGrowthTableJohnsonECCC,
                      modelVersion = "Johnson",
                      
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

  popGrowthPars <- demographicCoefficients(
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
  nr=10000
  cs = compositionBiasCorrection(w=modList$cowMult,q=runif(nr,modList$qMin,modList$qMax),
                                u=runif(nr,modList$uMin,modList$uMax),
                                z=runif(nr,modList$zMin,modList$zMax),approx=T)
  bias.Prior1 = cs$mu
  bias.Prior2 = cs$sig2^0.5
      
  if (returnValues) {
    betaPriors <- list(
      l.R.Prior1 = rPriorCoefs$Intercept,
      l.R.Prior2 = rPriorStdErrs$Intercept * modList$rIntSEMod,
      beta.Rec.anthro.Prior1 = rPriorCoefs$Anthro,
      beta.Rec.anthro.Prior2 = rPriorStdErrs$Anthro * modList$rAnthroSlopeSEMod,
      beta.Rec.fire.Prior1 = rPriorCoefs$fire_excl_anthro,
      beta.Rec.fire.Prior2 = rPriorStdErrs$fire_excl_anthro,
      sig.R.Prior1 = modList$rInterannualVar,
      sig.R.Prior2 = modList$rInterannualVarSE,
      l.Saf.Prior1 = sPriorCoefs$Intercept,
      l.Saf.Prior2 = sPriorStdErrs$Intercept * modList$sIntSEMod,
      beta.Saf.Prior1 = sPriorCoefs$Anthro,
      beta.Saf.Prior2 = sPriorStdErrs$Anthro * modList$sAnthroSlopeSEMod,
      sig.Saf.Prior1 = modList$sInterannualVar,
      sig.Saf.Prior2 = modList$sInterannualVarSE,
      bias.Prior1 = bias.Prior1,
      bias.Prior2 = bias.Prior2
    )

    # replace NULL values with 0 for when anthro or fire is not included
    betaPriors <- lapply(betaPriors, function(x) {
      if (is.null(x) || length(x) == 0) {
        1e-10
      } else {
        x
      }
    })
  } else {
    betaPriors <- list(
      l.R.Prior1 = rPriorCoefs$Intercept,
      l.R.Prior2 = paste0(round(rPriorStdErrs$Intercept, 4), "*", modList$rIntSEMod),
      beta.Rec.anthro.Prior1 = rPriorCoefs$Anthro,
      beta.Rec.anthro.Prior2 = paste0(round(rPriorStdErrs$Anthro, 4), "*",
                                      modList$rAnthroSlopeSEMod),
      beta.Rec.fire.Prior1 = rPriorCoefs$fire_excl_anthro,
      beta.Rec.fire.Prior2 = rPriorStdErrs$fire_excl_anthro,
      sig.R.Prior1 = modList$rInterannualVar,
      sig.R.Prior2 = modList$rInterannualVarSE,
      l.Saf.Prior1 = sPriorCoefs$Intercept,
      l.Saf.Prior2 = paste0(round(sPriorStdErrs$Intercept, 4), "*", modList$sIntSEMod),
      beta.Saf.Prior1 = sPriorCoefs$Anthro,
      beta.Saf.Prior2 = paste0(round(sPriorStdErrs$Anthro, 4), "*", modList$sAnthroSlopeSEMod),
      sig.Saf.Prior1 = modList$sInterannualVar,
      sig.Saf.Prior2 = modList$sInterannualVarSE,
      bias.Prior1 = bias.Prior1,
      bias.Prior2 = bias.Prior2
    )
  }
  return(betaPriors)
}

simCovariates <- function(initAnthro, initFire, numYears, anthroSlope,
                          anthroSlopeFuture, futureStep, fireSlope = 0) {
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

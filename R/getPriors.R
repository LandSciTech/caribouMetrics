#' Get model priors from national models
#'
#'
#'
#' @param modifiers a named list of modifiers to use to change the priors. If a
#'   modifier is supplied here the corresponding argument below is ignored.
#' @param returnValues logical. Default is TRUE. If FALSE returns strings for
#'   some values showing the initial values and the modifier ie "0.9 * 1.05"
#' @param bse anthropogenic disturbance slope survival uncertainty multiplier. 1
#'   - 10
#' @param bre anthropogenic disturbance slope recruitment uncertainty
#'   multiplier. 1 - 10
#' @param lse survival intercept uncertainty multiplier. 1 - 10
#' @param lre recruitment intercept uncertainty multiplier. 1 - 10
#' @param sse interannual coefficient of variation for survival. 0-1. See
#'   [popGrowthJohnson()] and functions therein for details
#' @param ssv uncertainty about interannual variation in survival. 0-1
#' @param sre interannual coefficient of variation for recruitment. 0-1
#' @param srv uncertainty about interannual variation in recruitment. 0-1.
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
#' * sig.Saf.Prior2: Uncertainty about interannual variation in adult female survival
#'
#' @examples
#' getPriors()
#' @export
getPriors <- function(modifiers = NULL,
                      survivalModelNumber = "M1",
                      recruitmentModelNumber = "M4",
                      bre = 4,
                      bse = 3,
                      lse = 5,
                      sse = 0.08696 * 0.4,
                      ssv = 0.03,
                      lre = 3,
                      sre = 0.46 * 0.5,
                      srv = 0.22,
                      populationGrowthTable = caribouMetrics::popGrowthTableJohnsonECCC,
                      modVer = "Johnson",
                      returnValues = TRUE) {
  # modifiers=cs

  expectMods <- c(as.list(environment()))
  expectMods$modifiers <- NULL
  if (is.null(modifiers)) {
    modifiers <- expectMods
  } else {
    # keep all values in modifiers and add any that are missing using values in
    # expectMods
    modifiers <- c(modifiers, expectMods[which(!names(expectMods) %in% names(modifiers))])
  }

  popGrowthPars <- demographicCoefficients(
    2,
    modelVersion = modVer,
    survivalModelNumber = modifiers$survivalModelNumber,
    recruitmentModelNumber = modifiers$recruitmentModelNumber,
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

  if (returnValues) {
    betaPriors <- list(
      l.R.Prior1 = rPriorCoefs$Intercept,
      l.R.Prior2 = rPriorStdErrs$Intercept * modifiers$lre,
      beta.Rec.anthro.Prior1 = rPriorCoefs$Anthro,
      beta.Rec.anthro.Prior2 = rPriorStdErrs$Anthro * modifiers$bre,
      beta.Rec.fire.Prior1 = rPriorCoefs$fire_excl_anthro,
      beta.Rec.fire.Prior2 = rPriorStdErrs$fire_excl_anthro,
      sig.R.Prior1 = modifiers$sre,
      sig.R.Prior2 = modifiers$srv,
      l.Saf.Prior1 = sPriorCoefs$Intercept,
      l.Saf.Prior2 = sPriorStdErrs$Intercept * modifiers$lse,
      beta.Saf.Prior1 = sPriorCoefs$Anthro,
      beta.Saf.Prior2 = sPriorStdErrs$Anthro * modifiers$bse,
      sig.Saf.Prior1 = modifiers$sse,
      sig.Saf.Prior2 = modifiers$ssv
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
      l.R.Prior2 = paste0(round(rPriorStdErrs$Intercept, 4), "*", modifiers$lre),
      beta.Rec.anthro.Prior1 = rPriorCoefs$Anthro,
      beta.Rec.anthro.Prior2 = paste0(round(rPriorStdErrs$Anthro, 4), "*",
                                      modifiers$bre),
      beta.Rec.fire.Prior1 = rPriorCoefs$fire_excl_anthro,
      beta.Rec.fire.Prior2 = rPriorStdErrs$fire_excl_anthro,
      sig.R.Prior1 = modifiers$sre,
      sig.R.Prior2 = modifiers$srv,
      l.Saf.Prior1 = sPriorCoefs$Intercept,
      l.Saf.Prior2 = paste0(round(sPriorStdErrs$Intercept, 4), "*", modifiers$lse),
      beta.Saf.Prior1 = sPriorCoefs$Anthro,
      beta.Saf.Prior2 = paste0(round(sPriorStdErrs$Anthro, 4), "*", modifiers$bse),
      sig.Saf.Prior1 = modifiers$sse,
      sig.Saf.Prior2 = modifiers$ssv
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

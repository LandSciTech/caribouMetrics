# this creates an environment where we can store objects that will be available
# to multiple functions/multiple function calls. Does not persist across
# sessions but it only take ~ 20s so once per session is probably ok.
# See explanation here: https://r-pkgs.org/data.html#sec-data-state
cacheEnv <- new.env()

if(file.exists("inst/extdata/simsNationalRadjusted.rds")){
  simsNationalRadjusted <- readRDS( "inst/extdata/simsNationalRadjusted.rds")
  simsNationalRunadjusted <- readRDS( "inst/extdata/simsNationalRunadjusted.rds")

  assign("simsNationalRadjusted", simsNationalRadjusted, envir = cacheEnv)
  assign("simsNationalRunadjusted", simsNationalRunadjusted, envir = cacheEnv)
}

#' Get a set of simulation results from the national model
#'
#'
#' @param Anthro,fire_excl_anthro numeric. A vector of numbers between 0 and 100
#'   representing the percentage of the landscape covered by anthropogenic
#'   disturbance buffered by 500 m, and the percentage covered by fire that does
#'   not overlap anthropogenic disturbance. The two vectors will be combined
#'   with `expand.grid()` to give the set of scenarios simulated.
#' @param forceUpdate logical. If the default inputs are used the result is
#'   cached. Set `forceUpdate` to TRUE to ensure the simulations are re-run.
#' @inheritParams demographicCoefficients
#' @inheritParams popGrowthJohnson
#' @param N0 initial population size
#'
#' @return a list with two elements:
#'  * summary: a tibble with a summary of parameter values for each scenario.
#'    Column names are Anthro, Mean, lower, upper, Parameter.
#'  * samples: a tibble with parameter values for each scenario and replicate
#'    4 rows per replicate \* scenario. Column names are Anthro, Parameter and Value
#' @export
#'
#' @examples
#' getSimsNational()
getSimsNational <- function(replicates = 1000, N0 = 1000, Anthro = seq(0, 100, by = 1),
                            fire_excl_anthro = 0, useQuantiles  = NULL,
                            populationGrowthTable  = NULL, adjustR = F, forceUpdate = F) {
  # replicates=1000;N0=1000;Anthro=seq(0,100,by=1);fire_excl_anthro=0;
  # useQuantiles =NULL;adjustR=F;forceUpdate=F
  doSave <- FALSE

  # check that everything other than adjustR is default
  check <- as.list(match.call())
  check$adjustR <- NULL

  saveName <- ifelse(adjustR, "simsNationalRadjusted", "simsNationalRunadjusted")

  if (length(check) == 1) {
    if (exists(saveName, envir=cacheEnv)) {
      message("Using saved object")
      return(get(saveName, envir=cacheEnv))
    } else {
      message("Object will be saved for future use")
      doSave <- TRUE
    }
  }

  check$forceUpdate <- NULL

  if (forceUpdate & (length(check) == 1)) {
    message("Updating cached national simulations.")
    doSave <- T
  }
  covTableObs <- expand.grid(
    Anthro = Anthro,
    fire_excl_anthro = fire_excl_anthro
  )
  covTableObs$Total_dist <- covTableObs$Anthro + covTableObs$fire_excl_anthro

  if (is.null(populationGrowthTable )) {
    populationGrowthTable  <- caribouMetrics::popGrowthTableJohnsonECCC
  }
  if (is.null(useQuantiles )) {
    popGrowthPars <- demographicCoefficients(
      replicates,
      populationGrowthTable = populationGrowthTable
    )
    rateSamplesAll <- demographicRates(covTable = covTableObs,
                                       popGrowthPars = popGrowthPars,
                                       returnSample = T, useQuantiles = F)
  } else {
    popGrowthPars <- demographicCoefficients(
      replicates, useQuantiles = useQuantiles,
      populationGrowthTable = populationGrowthTable
    )
    rateSamplesAll <- demographicRates(covTable = covTableObs,
                                       popGrowthPars = popGrowthPars,
                                       returnSample = T)
  }
  pars <- merge(data.frame(N0 = N0), rateSamplesAll)
  pars <- cbind(pars, popGrowthJohnson(pars$N0, R_bar = pars$R_bar,
                                       S_bar = pars$S_bar, numSteps = 1,
                                       K = F, adjustR = adjustR))
  simSurvBig <- pars %>%
    select("Anthro", "S_t") %>%
    group_by(Anthro) %>%
    summarize(Mean = mean(S_t), lower = quantile(S_t, 0.025),
              upper = quantile(S_t, 0.975))
  simSurvBig$Parameter <- "Adult female survival"
  simRecBig <- pars %>%
    select("Anthro", "R_t") %>%
    group_by(Anthro) %>%
    summarize(Mean = mean(R_t), lower = quantile(R_t, 0.025),
              upper = quantile(R_t, 0.975))
  simRecBig$Parameter <- "Recruitment"
  simLamBig <- pars %>%
    select("Anthro", "lambda") %>%
    group_by(Anthro) %>%
    summarize(Mean = mean(lambda), lower = quantile(lambda, 0.025),
              upper = quantile(lambda, 0.975))
  simLamBig$Parameter <- "Population growth rate"
  simFpopBig <- pars %>%
    select("Anthro", "N") %>%
    group_by(Anthro) %>%
    summarize(Mean = mean(N), lower = quantile(N, 0.025),
              upper = quantile(N, 0.975))
  simFpopBig$Parameter <- "Female population size"
  simBig <- rbind(simSurvBig, simRecBig, simLamBig, simFpopBig)

  parsSelect <- subset(pars, select = c(Anthro, S_t, R_t, lambda, N))
  names(parsSelect) <- c("Anthro", "Adult female survival",
                         "Recruitment", "Population growth rate",
                         "Female population size")
  parsSelect <- parsSelect %>%
    tidyr::pivot_longer(!Anthro, names_to = "Parameter", values_to = "Value")

  simBig <- list(summary = simBig, samples = parsSelect)

  if (doSave) {
    assign(saveName, simBig, envir = cacheEnv)
  }
  return(simBig)
}

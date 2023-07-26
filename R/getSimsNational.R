# this creates an environment where we can store objects that will be available
# to multiple functions/multiple function calls. Does not persist across
# sessions but it only take ~ 20s so once per session is probably ok.
# See explanation here: https://r-pkgs.org/data.html#sec-data-state
cacheEnv <- new.env(parent = emptyenv())

# Not supposed to save files to user computer on CRAN so for users the cache is
# only preserved within a session but for dev I have added this "persistent
# cache" use savePersistentCache function to update/create it after having run
# getSimsNational
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
#' @inheritParams caribouPopGrowth
#' @param N0 initial population size
#' @param cPars optional. Parameters for calculating composition survey bias term.
#'
#' @return a list with two elements:
#'  * summary: a tibble with a summary of parameter values for each scenario.
#'    Column names are Anthro, Mean, lower, upper, Parameter.
#'  * samples: a tibble with parameter values for each scenario and replicate
#'    4 rows per replicate \* scenario. Column names are Anthro, Parameter and Value
#' 
#' @family demography
#' @export
#'
#' @examples
#' getSimsNational()
getSimsNational <- function(replicates = 1000, N0 = 1000, Anthro = seq(0, 100, by = 1),
                            fire_excl_anthro = 0, useQuantiles  = NULL,
                            populationGrowthTable  = NULL, adjustR = TRUE, cPars=getScenarioDefaults(), forceUpdate = F) {
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
      doSave <- TRUE
    }
  }

  check$forceUpdate <- NULL

  if (forceUpdate & (length(check) == 1)) {
    doSave <- TRUE
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
                                       returnSample = TRUE, useQuantiles = FALSE)
  } else {
    popGrowthPars <- demographicCoefficients(
      replicates, useQuantiles = useQuantiles,
      populationGrowthTable = populationGrowthTable
    )
    rateSamplesAll <- demographicRates(covTable = covTableObs,
                                       popGrowthPars = popGrowthPars,
                                       returnSample = T)
  }
  
  bc = unique(subset(rateSamplesAll,select=replicate));nr=nrow(bc)
  bc$c = compositionBiasCorrection(q=runif(nr,cPars$qMin,cPars$qMax),w=cPars$cowMult,u=runif(nr,cPars$uMin,cPars$uMax),
                                   z=runif(nr,cPars$zMin,cPars$zMax))
  rateSamplesAll$c = NULL; rateSamplesAll= merge(rateSamplesAll, bc)
  
  pars <- merge(data.frame(N0 = N0), rateSamplesAll)
  pars <- cbind(pars, caribouPopGrowth(pars$N0, R_bar = pars$R_bar,
                                       S_bar = pars$S_bar, numSteps = 1,
                                       K = FALSE, adjustR = adjustR, c=pars$c, progress = FALSE))
  simSurvBig <- pars %>%
    select("Anthro", "S_t") %>%
    group_by(.data$Anthro) %>%
    summarize(Mean = mean(.data$S_t), lower = quantile(.data$S_t, 0.025),
              upper = quantile(.data$S_t, 0.975))
  simSurvBig$Parameter <- "Adult female survival"
  simRecBig <- pars %>%
    select("Anthro", "R_t") %>%
    group_by(.data$Anthro) %>%
    summarize(Mean = mean(.data$R_t), lower = quantile(.data$R_t, 0.025),
              upper = quantile(.data$R_t, 0.975))
  simRecBig$Parameter <- "Recruitment"
  simXBig <- pars %>%
    select("Anthro", "X_t") %>%
    group_by(.data$Anthro) %>%
    summarize(Mean = mean(.data$X_t), lower = quantile(.data$X_t, 0.025),
              upper = quantile(.data$X_t, 0.975))
  simXBig$Parameter <- "Adjusted recruitment"
  
  simLamBig <- pars %>%
    select("Anthro", "lambda") %>%
    group_by(.data$Anthro) %>%
    summarize(Mean = mean(.data$lambda), lower = quantile(.data$lambda, 0.025),
              upper = quantile(.data$lambda, 0.975))
  simLamBig$Parameter <- "Population growth rate"
  simFpopBig <- pars %>%
    select("Anthro", "N") %>%
    group_by(.data$Anthro) %>%
    summarize(Mean = mean(.data$N), lower = quantile(.data$N, 0.025),
              upper = quantile(.data$N, 0.975))
  simFpopBig$Parameter <- "Female population size"
  simBig <- rbind(simSurvBig, simRecBig, simXBig, simLamBig, simFpopBig)

  parsSelect <- subset(pars, select = c("Anthro", "S_t", "R_t","X_t", "lambda", "N"))
  names(parsSelect) <- c("Anthro", "Adult female survival",
                         "Recruitment","Adjusted recruitment", "Population growth rate",
                         "Female population size")
  parsSelect <- parsSelect %>%
    tidyr::pivot_longer(!.data$Anthro, names_to = "Parameter", values_to = "Value")

  simBig <- list(summary = simBig, samples = parsSelect)

  if (doSave) {
    message("Updating cached national simulations.")
    assign(saveName, simBig, envir = cacheEnv)
  }
  return(simBig)
}

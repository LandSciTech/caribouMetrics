test_that("testScript still works", {
  eParsIn <- list()
  eParsIn$collarOnTime <- 4
  eParsIn$collarOffTime <- 4
  eParsIn$collarNumYears <- 4

  scns <- expand.grid(
    obsYears = 8, collarCount = 30, cowMult = 2, collarInterval = 2,
    iAnthro = 0,
    tA = 0, obsAnthroSlope = 0, projAnthroSlope = 0, sQuantile = 0.960908218594268,
    rQuantile = 0.744425233039074, N0 = 1000
  )

  ##########
  # Get full set of sims for comparison
  simBig <- suppressWarnings(getSimsInitial(cPars = scns)) # If called with default parameters, use saved object to speed things up.

  ###############
  # Step 1: confirm appropriate prior variability in survival intercept using minimal (2) observed data points & 0 fire/anthro covariates. Controlled by priors on l.Saf, phi and sig.Saf.
  #################
  # source("CaribouDemoFns.R")
  # eParsIn$collarNumYears=1

  scResults <- suppressWarnings(runScnSet(scns,  simBig, eParsIn,
    niters = 100
  ))

  expect_s3_class(scResults$rr.summary.all, "data.frame")

  if (interactive()) {
    print(plotRes(scResults, "Population growth rate",lowBound = 0, highBound = 1.5))
    print(plotRes(scResults, "Recruitment"))
    print(plotRes(scResults, "Adult female survival"))
  }
})

test_that("bboutools scnenario with no disturbance and no additional monitoring ", {
  mod_flc <- here::here("results/test_mod_realc.rds")
  if (file.exists(mod_flc)) {
    mod_realc <- readRDS(mod_flc)
  } else {
    mod_realc <- bbouMakeSummaryTable(bboudata::bbousurv_a %>% filter(Year > 2010),
      bboudata::bbourecruit_a %>% filter(Year > 2010),
      N0 = NA, return_mcmc = T, niters = 3000
    )
    if(dir.exists(dirname(mod_flc))){
      saveRDS(mod_realc, mod_flc)
    }
  }

  simBig <- getSimsInitial(mod_realc)

  ###############
  # Example scenario - no disturbance and no additional monitoring
  scns <- data.frame(obsAnthroSlope = NA, projAnthroSlope = NA)
  scns$obsYears <- max(simBig$recruit_data$Year[!is.na(simBig$recruit_data$Calves)]) - min(simBig$recruit_data$Year) + 1
  scns$startYear <- min(simBig$recruit_data$Year)
  scns$projYears <- max(simBig$summary$Year) - scns$obsYears - scns$startYear
  scns$collarCount <- 0

  # devtools::load_all(path = "../caribouMetrics/")
  posteriorResult <- runScnSet(scns, simBig, niters = 3000)
  posteriorResult$obs.all <- NULL
  recPosterior <- plotRes(posteriorResult, "Recruitment")

  if (interactive()) {
    recPosterior
    # expect bands to match
  }

  # compare intervals
  recPosterior$data %>%
    select(-grp) %>%
    pivot_longer(c(Mean, lower, upper)) %>%
    pivot_wider(names_from = Type, values_from = value) %>%
    mutate(diff = abs(Bayesian - initial)) %>%
    pull(diff) %>%
    max() %>%
    # less than 1% absolute difference
    expect_lt(0.01)

  survPosterior <- plotRes(posteriorResult, "Adult female survival")

  if (interactive()) {
    survPosterior
    # expect bands to match
  }

  # compare intervals
  survPosterior$data %>%
    select(-grp) %>%
    pivot_longer(c(Mean, lower, upper)) %>%
    pivot_wider(names_from = Type, values_from = value) %>%
    mutate(diff = abs(Bayesian - initial)) %>%
    pull(diff) %>%
    max() %>%
    # less than 2% absolute difference
    expect_lt(0.02)
})

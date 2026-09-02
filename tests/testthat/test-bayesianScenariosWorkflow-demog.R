test_that("testScript still works", {
  eParsIn <- list()
  eParsIn$collarOnTime <- 4
  eParsIn$collarOffTime <- 4
  eParsIn$collarNumYears <- 4

  scns <- expand.grid(
    obsYears = c(8, 20), collarCount = 30, cowMult = 2, collarInterval = 2,
    iAnthro = 0,
    tA = 0, obsAnthroSlope = 0, projAnthroSlope = 0, sQuantile = 0.960908218594268,
    rQuantile = 0.744425233039074, N0 = 1000,N.sd=0.2
  )

  ##########
  # Get full set of sims for comparison
  simBig <- suppressWarnings(trajectoriesFromNational(cPars = scns)) # If called with default parameters, use saved object to speed things up.

  ###############
  # Step 1: confirm appropriate prior variability in survival intercept using minimal (2) observed data points & 0 fire/anthro covariates. Controlled by priors on l.Saf, phi and sig.Saf.
  #################
  # source("CaribouDemoFns.R")
  # eParsIn$collarNumYears=1

  scResults <- suppressWarnings(bayesianScenariosWorkflow(scns, simBig, eParsIn,
    niters = 100, printProgress = TRUE
  ))

  expect_s3_class(scResults$rr.summary.all, "data.frame")

  expect_warning(
    plotCompareTrajectories(scResults, "Recruitment"), 
    "duplicate"
  )
  expect_no_warning(
    plotCompareTrajectories(scResults, "Recruitment", facetVars = "ID")
  )
  
  if (interactive()) {
    print(plotCompareTrajectories(scResults, "Population growth rate", facetVars = "ID",
                                  lowBound = 0, highBound = 1.5))
    print(plotCompareTrajectories(scResults, "Recruitment", facetVars = "ID"))
    print(plotCompareTrajectories(scResults, "Adult female survival", facetVars = "ID"))
  }
  
  # test what happens if samples are returned from trajectoriesFromNational
  simBig2 <- trajectoriesFromNational(cPars = scns, returnSamples = TRUE)
  
  # # If scn table sets trajectory related parameters and simInitial has samples 
  # #   warn that simInitial$samples will be used.
   scResults2 <- expect_warning(
     bayesianScenariosWorkflow(scns, simBig2, eParsIn,
                               niters = 100, printProgress = TRUE))
  # 
   if(interactive()){
     plotCompareTrajectories(scResults2, "Population growth rate",
                             lowBound = 0, highBound = 1.5,facetVars="Replicate")
   }
})

test_that("Model and input trajectory match", {
  mod_flc <- here::here("results/test_bbou_nodist_nomonitor.rds")
  if (file.exists(mod_flc)) {
    mod_realc <- readRDS(mod_flc)
  } else {
    mod_realc <- estimateBayesianRates(
      bboudata::bbousurv_a %>% getCaribouYear() %>% filter(CaribouYear %>% between(2010, 2015)),
      bboudata::bbourecruit_a %>% getCaribouYear() %>% filter(CaribouYear %>% between(2010, 2015)),
      N0 = NA, return_mcmc = T, niters = 3000
    )
    if(dir.exists(dirname(mod_flc))){
      saveRDS(mod_realc, mod_flc)
    }
  }

  simBig <- trajectoriesFromBayesian(mod_realc)

  ###############
  # Example scenario - no disturbance and no additional monitoring
  scns <- data.frame(obsAnthroSlope = NA, projAnthroSlope = NA)
  inDat <- getCaribouYear(simBig$recruit_data)
  scns$obsYears <- max(inDat$CaribouYear[!is.na(inDat$Calves)]) - min(inDat$CaribouYear) + 1 
  scns$startYear <- min(inDat$CaribouYear)
  scns$projYears <- max(inDat$CaribouYear) - scns$obsYears - scns$startYear + 1
  scns$collarCount <- 0

  # devtools::load_all(path = "../caribouMetrics/")
  posteriorResult <- bayesianScenariosWorkflow(scns, simBig, niters = 3000,returnSamples=T)

  #for(nn in union(names(simBig$surv_data),names(posteriorResult$out$result$surv_data))){
  #  print(nn)
  #  expect_identical(simBig$surv_data[[nn]],posteriorResult$out$result$surv_data[[nn]])
  #}
  expect_identical(simBig$surv_data,posteriorResult$out$result$surv_data)
  
  posteriorResult$obs.all <- NULL
  recPosterior <- plotCompareTrajectories(posteriorResult, "Recruitment")

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

  survPosterior <- plotCompareTrajectories(posteriorResult, "Adult female survival")

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
    # less than 3% absolute difference
    expect_lt(0.03)
})


test_that("results match expected", {
  # save to speed up tests
  # simBig <- suppressWarnings(trajectoriesFromNational(N0 = 3000, forceUpdate=T))
  # saveRDS(simBig, "tests/testthat/data/simBig3000.rds", version = 2)
  
  simBig <- readRDS( file.path(test_path(), "data/simBig3000.rds"))
  doScn <- function(nCollar = 2000, nobsYears = 10, collarOn = 4, collarOff = 4, 
                    iAnthro = 0, obsAnthroSlope = 0, projAnthroSlope = 0, 
                    sQuantile = 0.5,  rQuantile = 0.5, rSlopeMod = 1, sSlopeMod = 1){
    #nCollar = 2000; nobsYears = 10; collarOn = 1; collarOff = 12; 
    #iAnthro = 0; obsAnthroSlope = 0; projAnthroSlope = 0; 
    #sQuantile = 0.5;  rQuantile = 0.5; rSlopeMod = 1; sSlopeMod = 1; KSDists = FALSE
    eParsIn <- list()
    eParsIn$collarOnTime <- collarOn
    eParsIn$collarOffTime <- collarOff
    eParsIn$collarNumYears <- 5
    
    scns <- expand.grid(
      obsYears = nobsYears, collarCount = nCollar, cowMult = 6, collarInterval = 1,
      assessmentYrs = 1, iAnthro = iAnthro, rSlopeMod = rSlopeMod, sSlopeMod = sSlopeMod,
      obsAnthroSlope = obsAnthroSlope, projAnthroSlope = projAnthroSlope,
      sQuantile = sQuantile, rQuantile = rQuantile, N0 = 3000,interannualVar=c("list(R_CV=0.23,S_CV=0.087)")
    )
    scResults <- suppressWarnings(bayesianScenariosWorkflow(
      scns, simBig,  eParsIn,
      niters = 3000
    ))
  }
  
  doPlot <- function(scResults, var = "Recruitment", title = "",highBound = 1){
    if (interactive()) {
      return(plotCompareTrajectories(scResults, var,
                                     lowBound = 0,  highBound = highBound,facetVars = NULL
      )+
        ggplot2::ggtitle(title))
    }
  }
  
  # difference between observed and true simulated observations
  calcDif <- function(obs, var){
    obs %>%
      filter(!MetricTypeID %in% c("Anthro", "Fire_excl_anthro")) %>% 
      select(Year, Mean, Type, Metric) %>% 
      tidyr::pivot_wider(names_from = "Type", values_from = "Mean") %>% 
      filter(!is.na(observed)) %>% 
      mutate(dif = abs(true - observed)) %>%
      group_by(Metric) %>% 
      summarise(mean_dif = mean(dif))
  }
  
  # difference between modeled and true simulated observations
  calcDifMod <- function(mod, var){
    #mod = manyObs
    obs_true <- mod$obs.all %>% 
      filter(!MetricTypeID %in% c("Anthro", "Fire_excl_anthro", "c"),
             Type == "true") %>% 
      select(Year, Mean, Type, Metric) 
    mod_proj <- mod$rr.summary.all %>% 
      mutate(ci_width = upper - lower, .keep = "unused") %>% 
      select(Year, Mean, Metric, ci_width)
    
    comp <- inner_join(obs_true, mod_proj, by = c("Year", "Metric"),
                       suffix = c("_true", "_proj")) %>% 
      # female pop size is done differently so don't compare
      filter(Metric != "Female population size") %>% 
      mutate(dif = Mean_true - Mean_proj) %>% 
      group_by(Metric) %>% 
      summarise(mean_dif = mean(abs(dif)),
                ci_width = mean(ci_width))
    comp
  }
  
  # difference between initial model and Bayesian model
  calcDifNat <- function(mod, min_year = 0){
    mod$rr.summary.all %>% select(Metric, Mean, Year) %>% 
      filter(Metric != "c") %>% 
      inner_join(mod$sim.all %>% select(Metric, Mean, Year),
                 by = c("Metric", "Year"), suffix = c("_PM", "_nat")) %>% 
      mutate(dif = Mean_PM - Mean_nat) %>% 
      filter(Year >= min_year) %>% 
      group_by(Metric) %>% 
      summarise(mean_dif = mean(dif)) %>% 
      # nat model does not do female adult pop in the same way
      filter(Metric != "Female population size")
  }
  
  # A pop with quantile >> 0.5 will be above the initial model projection
  highQ <- doScn(rQuantile = 0.95, sQuantile = 0.95)
  doPlot(highQ, "Adult female survival")
  
  difHighQ <- calcDifNat(highQ)
  
  expect_true(all(difHighQ$mean_dif > 0))
  
  # a pop that is less sensitive to anthro dist ie r/sSlopeMod < 1 will show a
  # line that diverges from the initial model. But only if there was some
  # disturbance in training data?
  lowSens <- doScn(rSlopeMod = 0.1, sSlopeMod = 0.1, iAnthro = 80, nobsYears = 20,
                   obsAnthroSlope = 1, projAnthroSlope = 1)
  doPlot(lowSens)
  doPlot(lowSens, "Adult female survival")
  
  difLowSens <- calcDifNat(lowSens, 2040)
  
  expect_true(all(difLowSens %>% pull(mean_dif) > 0))
  
  # same but no anthro in training data
  lowSensNtrain <- doScn(rSlopeMod = 0.1, sSlopeMod = 0.1, iAnthro = 0, nobsYears = 20,
                         obsAnthroSlope = 0, projAnthroSlope = 10)
  
  doPlot(lowSensNtrain)
  doPlot(lowSensNtrain, "Adult female survival")
  doPlot(lowSensNtrain, "Population growth rate",highBound=1.5)
  difLowSensNtrain <- calcDifNat(lowSensNtrain, min_year = 2040)
  
  # When no disturbance in training data obs don't tell the model that the
  # population is not sensitive to disturbance so it follows the national priors
  # but when there is disturbance in the training data the projections are very
  # different from the national model
  expect_true(all(difLowSensNtrain$mean_dif < difLowSens$mean_dif))
})
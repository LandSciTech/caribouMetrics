
test_that("Runs with defaults", {
  # reduce some defaults to make fast
  # note that the default csv does not match the default startYear
  expect_s3_class(caribouBayesianIPM(startYear = 2009, Nchains = 1, Niter = 100, Nburn = 10,
                       Nthin = 2)$result,
            "rjags")
})

test_that("input tables are as expected",{
  # default input data
  survDataIn <- system.file("extdata/simSurvData.csv",
                            package = "caribouMetrics") %>%
    read.csv()
  ageRatioIn <- system.file("extdata/simAgeRatio.csv",
                                 package = "caribouMetrics")%>%
    read.csv()
  disturbanceIn <- system.file("extdata/simDisturbance.csv",
                               package = "caribouMetrics")%>%
    read.csv()

  # default start year is outside range of data but still runs
  res1 <- expect_warning(caribouBayesianIPM(Nchains = 1, Niter = 100, Nburn = 10, Nthin = 2))

  expect_true(is.na(res1$survInput$surv[1]))
  expect_s3_class(res1$result, "rjags")

  # end year is outside range of data but still runs
  res2 <- expect_warning(caribouBayesianIPM(startYear = 2009, endYear = 2050, Nchains = 1, Niter = 100, Nburn = 10,
             Nthin = 2))

  expect_true(is.na(last(res2$survInput$surv)))
  expect_s3_class(res2$result, "rjags")

  # ageRatio is outside year range warning but still runs
  res3 <- expect_warning(caribouBayesianIPM(ageRatio = mutate(ageRatioIn, Year = Year - 30),
             startYear = 2009, endYear = 2040, Nchains = 1, Niter = 100, Nburn = 10,
             Nthin = 2), "composition")

  expect_s3_class(res3$result, "rjags")

  # all survival data is outside year range
  expect_error(
    caribouBayesianIPM(survData = mutate(survDataIn, Year = Year - 30),
               startYear = 2009, endYear = 2040, Nchains = 1, Niter = 100, Nburn = 10,
               Nthin = 2),
    "None of the survival data")

  # survival should go to curYear but disturbance data should go to endYear
  # no warnings when survData is shorter
  expect_s3_class(
    caribouBayesianIPM(survData = survDataIn,
               startYear = 2009, endYear = 2040, Nchains = 1, Niter = 100, Nburn = 10,
               Nthin = 2)$result,
    "rjags"
  )

  # warning when survData is missing a year in the middle
  expect_warning(
    expect_s3_class(
      caribouBayesianIPM(survData = filter(survDataIn, Year != 2010),
                 startYear = 2009, endYear = 2040, Nchains = 1, Niter = 100, Nburn = 10,
                 Nthin = 2)$result,
      "rjags"
    ),
    "consecutive years")

  # all disturbance data is outside year range
  expect_warning(
    expect_error(
      caribouBayesianIPM(disturbance = mutate(disturbanceIn, Year = Year - 50),
                 startYear = 2009, endYear = 2040, Nchains = 1, Niter = 100, Nburn = 10,
                 Nthin = 2),
      "None of the disturbance data")
  )

  # wrong column names
  expect_error(
    caribouBayesianIPM(disturbance = rename(disturbanceIn, year = Year),
               startYear = 2009, Nchains = 1, Niter = 100, Nburn = 10,
               Nthin = 2),
    "missing expected columns")

  expect_error(
    caribouBayesianIPM(ageRatio =  rename(ageRatioIn, cls = Class),
               startYear = 2009, Nchains = 1, Niter = 100, Nburn = 10,
               Nthin = 2),
    "missing expected columns")

  expect_error(
    caribouBayesianIPM(survData = rename(survDataIn, events = enter),
               startYear = 2009, Nchains = 1, Niter = 100, Nburn = 10,
               Nthin = 2),
    "missing expected columns")

  # check haven't added need for Total_dist
  expect_s3_class(
    caribouBayesianIPM(disturbance = select(disturbanceIn, -Total_dist),
               startYear = 2009, Nchains = 1, Niter = 100, Nburn = 10,
               Nthin = 2)$result,
    "rjags")
})

test_that("survAnalysisMethod works", {
  expect_message(out1 <- caribouBayesianIPM(startYear = 2009, Nchains = 1, Niter = 100, Nburn = 10,
                            Nthin = 2),
                 "using Kaplan-Meier survival model")
  expect_s3_class(out1$result, "rjags")

  expect_message(out2 <- caribouBayesianIPM(startYear = 2009, Nchains = 1, Niter = 100, Nburn = 10,
                            Nthin = 2, survAnalysisMethod = "Exponential"),
                 "expanding survival record")

  expect_s3_class(out2$result, "rjags")

  expect_message(out3 <- caribouBayesianIPM(startYear = 2009, Nchains = 1, Niter = 100, Nburn = 10,
                            Nthin = 2, survAnalysisMethod = "other"),
                 "expanding survival record")

  expect_s3_class(out3$result, "rjags")
})

test_that("works when only 1 collared animal",{
  cowCounts <- data.frame(
    Year = 2012:2023,
    Count = 10,
    Class = "cow"
  )

  freqStartsByYear <- data.frame(
    Year = 2012:2023,
    numStarts = 1
  )

  scns <- getScenarioDefaults(obsYears = 12, collarCount = 1)

  oo <- simulateObservations(scns, cowCounts = cowCounts,
                             freqStartsByYear = freqStartsByYear)

  # ensure some deaths so it uses KM still
  oo$simSurvObs$event[3] <- 1

  expect_warning(
    out <- caribouBayesianIPM(
      survData = oo$simSurvObs, ageRatio = oo$ageRatioOut,
      disturbance = oo$simDisturbance,
      startYear = 2012, endYear = 2043,
      Nchains = 1, Niter = 100, Nburn = 10,
      Nthin = 2),
    "low sample size")

  expect_s3_class(out$result, "rjags")

  # FIXED: error seems to happen when there are few collars but enough data for KM.
  # This was happening in some of the paper simulations too.
  # Error in `$<-.data.frame`(`*tmp*`, tau, value = c(18.0902177774918, 17.7911133740675,  :
  #                                                     replacement has 7 rows, data has 4

})

test_that("results match expected", {
  simBig <- suppressWarnings(getSimsNational(N0 = 3000))
  
  doScn <- function(nCollar = 2000, nobsYears = 10, collarOn = 1, collarOff = 12, 
                    iAnthro = 0, obsAnthroSlope = 0, projAnthroSlope = 0, 
                    sQuantile = 0.5,  rQuantile = 0.5, rSlopeMod = 1, sSlopeMod = 1){
    eParsIn <- list()
    eParsIn$cowCounts <- data.frame(
      Year = 1981:2023,
      Count = nCollar,
      Class = "cow"
    )
    eParsIn$freqStartsByYear <- data.frame(
      Year = 1981:2023,
      numStarts = nCollar
    )
    eParsIn$collarOnTime <- collarOn
    eParsIn$collarOffTime <- collarOff
    eParsIn$collarNumYears <- 5
    
    scns <- expand.grid(
      obsYears = nobsYears, collarCount = nCollar, cowMult = 2, collarInterval = 1,
      assessmentYrs = 1, iAnthro = iAnthro, rSlopeMod = rSlopeMod, sSlopeMod = sSlopeMod,
      tA = 0, obsAnthroSlope = obsAnthroSlope, projAnthroSlope = projAnthroSlope,
      sQuantile = sQuantile, rQuantile = rQuantile, N0 = 3000
    )
    scResults <- suppressWarnings(runScnSet(scns, eParsIn, simBig, getKSDists = F))
  }
  
  doPlot <- function(scResults, var = "Recruitment"){
    if (interactive()) {
      print(plotRes(scResults$rr.summary.all, var,
                    obs = scResults$obs.all,
                    lowBound = 0, simRange = scResults$sim.all, facetVars = NULL
      ))
    }
  }
  
  # difference between observed and true simulated observations
  calcDif <- function(obs, var){
    obs %>%
      select(Year, Mean, type, parameter) %>% 
      tidyr::pivot_wider(names_from = "type", values_from = "Mean") %>% 
      filter(!is.na(observed)) %>% 
      mutate(dif = abs(true - observed)) %>%
      group_by(parameter) %>% 
      summarise(mean_dif = mean(dif))
  }
  
  # difference between observed and true simulated observations
  calcDifMod <- function(mod, var){
    obs_true <- mod$obs.all %>% filter(type == "true") %>% 
      select(Year, Mean, type, parameter) 
    mod_proj <- mod$rr.summary.all %>% 
      select(Year, Mean, Parameter)
    
    left_join(obs_true, mod_proj, by = c("Year", parameter = "Parameter"),
              suffix = c("_true", "_proj")) %>% 
      mutate(dif = Mean_true - Mean_proj) %>% 
      group_by(parameter) %>% 
      summarise(mean_dif = mean(abs(dif)))
  }
  
  # difference between national model and IPM
  calcDifNat <- function(mod, min_year = 0){
    mod$rr.summary.all %>% select(Parameter, Mean, Year) %>% 
      right_join(mod$sim.all %>% select(parameter, Mean, Year),
                 by = c(Parameter = "parameter", "Year"), suffix = c("_IPM", "_nat")) %>% 
      mutate(dif = Mean_IPM - Mean_nat) %>% 
      filter(Year >= min_year) %>% 
      group_by(Parameter) %>% 
      summarise(mean_dif = mean(dif)) %>% 
      # nat model does not do female adult pop in the same way
      filter(Parameter != "Female population size")
  }
  
  # expectations that make sense based on what we know and ensure results
  # are similar to expected
  
  # when we have a lot of collars the distance between observations and "true"
  # pop is smaller than when we have few.
  manyObs <- doScn()
  doPlot(manyObs)
  
  fewCollarObs <- doScn(nCollar = 30)
  doPlot(fewCollarObs)
  
  difMany <- calcDif(manyObs$obs.all)
  difFew <- calcDif(fewCollarObs$obs.all)
  expect_true(all(difFew$mean_dif > difMany$mean_dif))
  
  # model predictions should also be closer to true with more collar data
  modDifMany <- calcDifMod(manyObs)
  modDifFew <- calcDifMod(fewCollarObs)

  expect_true(all(modDifFew$mean_dif > modDifMany$mean_dif))
  
  # difference between modeled and true does not change much if collar on/off times
  # are different but still a year apart
  collOff3On4 <- doScn(collarOn = 4, collarOff = 3)
  doPlot(collOff3On4)

  modDifcollOff3On4 <- calcDifMod(collOff3On4)
  
  # not more than 5 times difference
  expect_true(all((modDifcollOff3On4$mean_dif - 
                     modDifMany$mean_dif)/modDifMany$mean_dif < 5))

  # When there are gaps in survival data bc collars on/off at different times
  collOff12On4 <- doScn(collarOn = 4, collarOff = 12)
  doPlot(collOff12On4)
  doPlot(collOff12On4, "Adult female survival")
  

  collOff6On4 <- doScn(collarOn = 4, collarOff = 6)
  doPlot(collOff6On4)
  
  obsDifcollOff12On4 <- calcDif(collOff12On4$obs.all)
  
  modDifcollOff12On4 <- calcDifMod(collOff12On4)
  
  # TODO: what is expectation here?
  # expect survival to be the same with same parameters 
  # might ignore first year of collar data if year 1 starts later than jan
  
  # K-M gets weird if there are any months in the middle without data JH avoided
  # this by dropping 1st year
  
  # Can look at K-M directly with use survInput from the model and then pass to
  # getKMestimates should be insensitive to when things fall off or time between
  # on and off
  obs <- simulateObservations(paramTable = getScenarioDefaults(obsYears = 12),
    cowCounts = data.frame(Year = 2012:2023, Count = 1000, Class = "cow"),
    freqStartsByYear = data.frame(Year = 2012:2023, numStarts = 100)
  )
  
  mod <- caribouBayesianIPM(obs$simSurvObs, obs$ageRatioOut, obs$simDisturbance,
                            startYear = 2012, endYear = 2023)
  
  getKMSurvivalEstimates(mod$survInput)
  
  #TODO finish this
  
  # A pop with quantile >> 0.5 will be above the national model projection
  highQ <- doScn(rQuantile = 0.95, sQuantile = 0.95)
  doPlot(highQ, "Adult female survival")
  
  difHighQ <- calcDifNat(highQ)
  
  expect_true(all(difHighQ$mean_dif > 0))
  
  # a pop that is less sensitive to anthro dist ie r/sSlopeMod < 1 will show a
  # line that diverges from the national model. But only if there was some
  # disturbance in training data?
  lowSens <- doScn(rSlopeMod = 0.2, sSlopeMod = 0.2, iAnthro = 5, nobsYears = 20,
                   obsAnthroSlope = 2, projAnthroSlope = 2)
  doPlot(lowSens)
  doPlot(lowSens, "Adult female survival")
  
  difLowSens <- calcDifNat(lowSens, 2023)
  
  expect_true(all(difLowSens$mean_dif > 0))
  
  # same but no anthro in training data
  lowSensNtrain <- doScn(rSlopeMod = 0.2, sSlopeMod = 0.2, iAnthro = 0, nobsYears = 20,
                         obsAnthroSlope = 0, projAnthroSlope = 4)
  
  doPlot(lowSensNtrain)
  doPlot(lowSensNtrain, "Adult female survival")
  doPlot(lowSensNtrain, "Population growth rate")
  difLowSensNtrain <- calcDifNat(lowSensNtrain, min_year = 2023)
  
  # expect differences to be small
  expect_true(all(difLowSensNtrain$mean_dif < difLowSens$mean_dif))
  
  # KS distances JH added to characterize deviation from national model bands
  # not just the mean. So should test what happens when there is no sample info
  # provided. distribution of means from national model vs IPM. Set standards
  # that when no obs the differences don't get much worse than they are now. See
  # the doc JH will send to show what we are looking for
  
  # How to scan for problems with convergence ie the Rhats from the JAGS model
})

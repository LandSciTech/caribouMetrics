
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
  
  cowCounts <- data.frame(
    Year = 2014:2023,
    Count = 2000,
    Class = "cow"
  )
  
  freqStartsByYear <- data.frame(
    Year = 2014:2023,
    numStarts = 2000
  )
  
  scns <- getScenarioDefaults(obsYears = 10, collarCount = 2000)
  
  oo <- simulateObservations(scns, cowCounts = cowCounts,
                             freqStartsByYear = freqStartsByYear)
  
  
  
  out <- caribouBayesianIPM(
    survData = oo$simSurvObs, ageRatio = oo$ageRatioOut,
    disturbance = oo$simDisturbance,
    startYear = scns$startYear, endYear = scns$curYear+scns$projYears,
    Nchains = 1, Niter = 100, Nburn = 10,
    Nthin = 2)
  
  # TODO: add expectations that make sense based on what we know and ensure obs
  # are similar to expected
})

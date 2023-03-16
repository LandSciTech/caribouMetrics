
test_that("Runs with defaults", {
  # reduce some defaults to make fast
  # note that the default csv does not match the default startYear
  expect_is(runRMModel(startYear = 2009, Nchains = 1, Niter = 100, Nburn = 10,
                       Nthin = 2),
            "list")
})

test_that("input tables are as expected",{
  # default input data
  survDataIn <- system.file("extdata/simSurvData.csv",
                            package = "BayesianCaribouDemographicProjection") %>%
    read.csv()
  ageRatio.herdIn <- system.file("extdata/simAgeRatio.csv",
                                 package = "BayesianCaribouDemographicProjection")%>%
    read.csv()
  disturbanceIn <- system.file("extdata/simDisturbance.csv",
                               package = "BayesianCaribouDemographicProjection")%>%
    read.csv()

  # default start year is outside range of data but still runs
  res1 <- expect_warning(runRMModel(Nchains = 1, Niter = 100, Nburn = 10, Nthin = 2))

  expect_true(is.na(res1$survInput$surv[1]))

  # end year is outside range of data but still runs
  res2 <- expect_warning(runRMModel(startYear = 2009, endYear = 2050, Nchains = 1, Niter = 100, Nburn = 10,
             Nthin = 2))

  expect_true(is.na(last(res2$survInput$surv)))

  # ageRatio.herd is outside year range warning but still runs
  expect_warning(runRMModel(ageRatio.herd = mutate(ageRatio.herdIn, Year = Year - 30),
             startYear = 2009, endYear = 2040, Nchains = 1, Niter = 100, Nburn = 10,
             Nthin = 2), "composition")

  # all survival data is outside year range
  expect_warning(
    expect_error(
      runRMModel(survData = mutate(survDataIn, Year = Year - 30),
                 startYear = 2009, endYear = 2040, Nchains = 1, Niter = 100, Nburn = 10,
                 Nthin = 2),
      "None of the survival data")
  )

  # all disturbance data is outside year range
  expect_warning(
    expect_error(
      runRMModel(disturbance = mutate(disturbanceIn, Year = Year - 50),
                 startYear = 2009, endYear = 2040, Nchains = 1, Niter = 100, Nburn = 10,
                 Nthin = 2),
      "None of the disturbance data")
  )

  # wrong column names
  expect_error(
    runRMModel(disturbance = rename(disturbanceIn, year = Year),
               startYear = 2009, Nchains = 1, Niter = 100, Nburn = 10,
               Nthin = 2),
    "missing expected columns")

  expect_error(
    runRMModel(ageRatio.herd =  rename(ageRatio.herdIn, cls = Class),
               startYear = 2009, Nchains = 1, Niter = 100, Nburn = 10,
               Nthin = 2),
    "missing expected columns")

  expect_error(
    runRMModel(survData = rename(survDataIn, events = enter),
               startYear = 2009, Nchains = 1, Niter = 100, Nburn = 10,
               Nthin = 2),
    "missing expected columns")

  # check haven't added need for Total_dist
  expect_is(
    runRMModel(disturbance = select(disturbanceIn, -Total_dist),
               startYear = 2009, Nchains = 1, Niter = 100, Nburn = 10,
               Nthin = 2),
    "list")
})

test_that("survAnalysisMethod works", {
  expect_message(runRMModel(startYear = 2009, Nchains = 1, Niter = 100, Nburn = 10,
                            Nthin = 2),
                 "using Kaplan-Meier survival model")
  # TODO fix this error
  expect_message(runRMModel(startYear = 2009, Nchains = 1, Niter = 100, Nburn = 10,
                            Nthin = 2, survAnalysisMethod = "Exponential"),
                 "using parametric exponential survival model")

  expect_message(runRMModel(startYear = 2009, Nchains = 1, Niter = 100, Nburn = 10,
                            Nthin = 2, survAnalysisMethod = "other"),
                 "expanding survival record")
})

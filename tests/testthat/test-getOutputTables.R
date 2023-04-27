test_that("works with defaults", {
  
  scns <- getScenarioDefaults(J = 0)
  simO <- simulateObservations(scns,
                               freqStartsByYear = data.frame(Year = 2014:2023,
                                                             numStarts = 20),
                               cowCounts = data.frame(Year = 2014:2023,
                                                      Count = 100,
                                                      Class = "cow"))
  
  out <- caribouBayesianIPM(survData = simO$simSurvObs, ageRatio.herd = simO$ageRatioOut,
                    disturbance = simO$simDisturbance,
                    startYear = 2014, Nchains = 1, Niter = 100, Nburn = 10,
                    Nthin = 2)
  
  # error when result has different startYear from argument
  expect_error(getOutputTables(out, startYear = 2009, endYear = 2023, oo = simO), 
               "different length")

  expect_type(getOutputTables(out, startYear = 2014, endYear = 2023, oo = simO,
                              simBig = getSimsNational(), getKSDists = FALSE), 
              "list")
})

test_that("works with defaults", {
  
  scns <- getScenarioDefaults(projYears = 0, obsYears = 10)
  simO <- simulateObservations(scns,
                               freqStartsByYear = data.frame(Year = 2014:2023,
                                                             numStarts = 20),
                               cowCounts = data.frame(Year = 2014:2023,
                                                      Count = 100,
                                                      Class = "cow"))
  
  out <- caribouBayesianIPM(survData = simO$simSurvObs, ageRatio = simO$ageRatioOut,
                    disturbance = simO$simDisturbance,
                    startYear = 2014, Nchains = 1, Niter = 100, Nburn = 10,
                    Nthin = 2)
  
  # error when result has different startYear from argument
  expect_error(getOutputTables(out, startYear = 2009, endYear = 2023, simObsList = simO), 
               "different length")

  expect_type(getOutputTables(out, startYear = 2014, endYear = 2023, simObsList = simO,
                              simNational = getSimsNational(), getKSDists = FALSE), 
              "list")
})

test_that("decimals in observed disturbance work", {
  scns <- getScenarioDefaults(projYears = 10, obsYears = 10, obsAnthroSlope = 1.5, projAnthroSlope = 5)
  simO <- simulateObservations(scns,
                               freqStartsByYear = data.frame(Year = 2014:2023,
                                                             numStarts = 20),
                               cowCounts = data.frame(Year = 2014:2023,
                                                      Count = 100,
                                                      Class = "cow"))
  
  out <- caribouBayesianIPM(survData = simO$simSurvObs, ageRatio = simO$ageRatioOut,
                            disturbance = simO$simDisturbance,
                            startYear = 2014, Nchains = 1, Niter = 100, Nburn = 10,
                            Nthin = 2)
  
  expect_message(getOutputTables(out, startYear = 2014, endYear = 2023, simObsList = simO, 
                  simNational = getSimsNational(), getKSDists = FALSE),
                 "recalculating")
             
})

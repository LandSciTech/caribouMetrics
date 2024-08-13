test_that("works with defaults", {
  
  scns <- getScenarioDefaults(projYears = 0, obsYears = 10, collarCount = 20,
                              cowMult = 3)
  simO <- simulateObservations(scns)
  
  out <- caribouBayesianPM(survData = simO$simSurvObs, ageRatio = simO$ageRatioOut,
                    disturbance = simO$simDisturbance,
                    Nchains = 1, Niter = 100, Nburn = 10,
                    Nthin = 2)
  
  # error when result has different startYear from argument
  expect_error(getOutputTables(out, startYear = 2009, endYear = 2023), 
               "different length")

  expect_type(getOutputTables(out, simNational = getSimsNational(),
                              getKSDists = FALSE), 
              "list")
})

test_that("decimals in observed disturbance work", {
  scns <- getScenarioDefaults(projYears = 10, obsYears = 10, 
                              obsAnthroSlope = 1.5, projAnthroSlope = 5,
                              collarCount = 20, cowMult = 3)
  simO <- simulateObservations(scns)
  
  out <- caribouBayesianPM(survData = simO$simSurvObs, ageRatio = simO$ageRatioOut,
                            disturbance = simO$simDisturbance,
                            Nchains = 1, Niter = 100, Nburn = 10,
                            Nthin = 2)
  
  expect_message(
    getOutputTables(
      out, exData = simO$exData, paramTable = simO$paramTable,
      simNational = getSimsNational(), getKSDists = FALSE),
    "recalculating")
             
})

test_that("works with out sim obs",{
  mod_real <- caribouBayesianPM(Niter = 100, Nburn = 10)
  
  
  mod_tbl <- getOutputTables(mod_real,
                             simNational = getSimsNational(),
                             getKSDists = FALSE)
  
  expect_type(mod_tbl, "list")
  
})

test_that("works with out simNational", {
  mod_real <- caribouBayesianPM(Niter = 100, Nburn = 10)
  
  
  mod_tbl <- getOutputTables(mod_real)
  
  expect_type(mod_tbl, "list")
})

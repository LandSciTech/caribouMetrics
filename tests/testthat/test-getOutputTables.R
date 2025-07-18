test_that("works with defaults", {
  
  scns <- getScenarioDefaults(projYears = 0, obsYears = 10, collarCount = 20,
                              cowMult = 3)
  trajs <- getSimsInitial()$samples
  simO <- simulateObservations(trajs, scns)
  
  out <- caribouBayesianPM(survData = simO$simSurvObs, recruitData = simO$simRecruitObs,
                    disturbance = simO$simDisturbance,
                    niters=100)
  
  # error when result has different startYear from argument
  expect_error(getOutputTables(out, startYear = 2009, endYear = 2023), 
               "different length")

  expect_type(getOutputTables(out, simInitial = getSimsInitial()), 
              "list")
})

test_that("decimals in observed disturbance work", {
  scns <- getScenarioDefaults(projYears = 10, obsYears = 10, 
                              obsAnthroSlope = 1.5, projAnthroSlope = 5,
                              collarCount = 20, cowMult = 3)
  simO <- simulateObservations(scns)
  
  out <- caribouBayesianPM(survData = simO$simSurvObs, recruitData = simO$simRecruitObs,
                            disturbance = simO$simDisturbance,
                            niters=100)
  
  expect_message(
    getOutputTables(
      out, exData = simO$exData, paramTable = simO$paramTable,
      simInitial = getSimsInitial()),
    "recalculating")
             
})

mod_real <- caribouBayesianPM(survData = bboudata::bbousurv_a %>% filter(Year > 2010), 
                              recruitData = bboudata::bbourecruit_a %>% filter(Year > 2010),
                              niters=1)

test_that("works with out sim obs",{
  mod_tbl <- getOutputTables(mod_real,
                             simInitial = getSimsInitial())
  
  expect_type(mod_tbl, "list")
  
})

test_that("works with out simInitial", {
  
  mod_tbl <- getOutputTables(mod_real)
  
  expect_type(mod_tbl, "list")
})

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
  # TODO: See Handle this case at simulateObservations.R#137
  # scns <- getScenarioDefaults(projYears = 10, obsYears = 10, 
  #                             obsAnthroSlope = 1.5, projAnthroSlope = 5,
  #                             collarCount = 20, cowMult = 3)
  # trajs <- getSimsInitial()$samples
  # simO <- simulateObservations(trajs, scns)
  # 
  # out <- caribouBayesianPM(survData = simO$simSurvObs, recruitData = simO$simRecruitObs,
  #                           disturbance = simO$simDisturbance,
  #                           niters=100)
  # 
  # expect_message(
  #   getOutputTables(
  #     out, exData = simO$exData, paramTable = simO$paramTable,
  #     simInitial = getSimsInitial()),
  #   "recalculating")
             
})

mod_real <- caribouBayesianPM(survData = bboudata::bbousurv_a %>% filter(Year > 2010), 
                              recruitData = bboudata::bbourecruit_a %>% filter(Year > 2010),
                              niters=1)

test_that("works with simInitial",{
  # Can't compare with out disturbance in the original model
  expect_error(getOutputTables(mod_real,
                               simInitial = getSimsInitial()), "Set disturbance")
  
  # Works when disturbance is specified
  disturbance <- unique(subset(bboudata::bbourecruit_a %>% filter(Year > 2010), 
                               select = Year))
  disturbance$Anthro <- 0
  disturbance$fire_excl_anthro <- 0
  mod_realb <- caribouBayesianPM(
    survData = bboudata::bbousurv_a %>% filter(Year > 2010),
    recruitData = bboudata::bbourecruit_a %>% filter(Year > 2010),
    disturbance = disturbance,
    niters = 10
  )
  
  mod_tbl <- getOutputTables(mod_realb,
                              simInitial = getSimsInitial())

  expect_type(mod_tbl, "list")
  
})

test_that("works with out simInitial", {
  
  mod_tbl <- getOutputTables(mod_real)
  
  expect_type(mod_tbl, "list")
})

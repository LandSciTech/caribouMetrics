test_that("real and simulated data work", {
  # Using observed survival, recruitment and disturbance data
  mod <- bayesianTrajectoryWorkflow(
    surv_data = bboudata::bbousurv_a %>% filter(Year > 2010),
    recruit_data = bboudata::bbourecruit_a %>% filter(Year > 2010),
    disturbance = NULL
  )
  str(mod, max.level = 2)
  
  expect_type(mod, "list")
  
  # Using simulated observation data
  scns <- getScenarioDefaults(projYears = 10, obsYears = 10,
                              obsAnthroSlope = 1, projAnthroSlope = 5,
                              collarCount = 20, cowMult = 5)
  
  simO <- simulateObservations(scns)
  
  out <- bayesianTrajectoryWorkflow(surv_data = simO$simSurvObs, recruit_data = simO$simRecruitObs,
                           disturbance = simO$simDisturbance,
                           startYear = 2014)
  
  expect_type(out, "list")
})

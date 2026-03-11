test_that("works with defaults", {
  #devtools::load_all()
  
  scns <- getScenarioDefaults(projYears = 0, obsYears = 10, collarCount = 0,
                              cowMult = 3)
  simIni <-trajectoriesFromNational() 
  simO <- simulateObservations(scns)
  
  out <- bayesianTrajectoryWorkflow(surv_data = simO$simSurvObs, 
                           recruit_data = simO$simRecruitObs,
                           disturbance = simO$simDisturbance,
                           niters=100)
  
  out_tbls <- compareTrajectories(out, simInitial = simIni, paramTable = scns)
  
  # all components are present
  purrr::map_lgl(out_tbls, \(x) nrow(x) > 0) %>% all() %>% 
    expect_true()
  
  #visual checks
  if(0){
    recPrior =  plotCompareTrajectories(out_tbls, "Recruitment")
    plot(recPrior)
    survPrior =  plotCompareTrajectories(out_tbls, "Adult female survival")
    plot(survPrior)
    lamPrior =  plotCompareTrajectories(out_tbls, "Population growth rate",lowBound=0.5,highBound=1.5)
    plot(lamPrior)
  }  
  #TO DO: confirm that out bands approximately match simIni bands - see visual checks.
})

#Note: see test-bayesianScenariosWorkflow.R for more output checks.

test_that("decimals in observed disturbance work", {
  # TODO: See Handle this case at simulateObservations.R#137
  # scns <- getScenarioDefaults(projYears = 10, obsYears = 10,
  #                             obsAnthroSlope = 1.5, projAnthroSlope = 5,
  #                             collarCount = 20, cowMult = 3)
  # simO <- simulateObservations(scns)
  # 
  # out <- bayesianTrajectoryWorkflow(surv_data = simO$simSurvObs, recruit_data = simO$simRecruitObs,
  #                           disturbance = simO$simDisturbance,
  #                           niters=100)
  # 
  # expect_message(
  #   compareTrajectories(
  #     out, exData = simO$exData, paramTable = simO$paramTable,
  #     simInitial = trajectoriesFromNational()),
  #   "recalculating")
             
})

# Use saved file because this takes a long time. 
mod_fl <- here::here("results/test_mod_real.rds")
if(file.exists(mod_fl)){
  mod_real <- readRDS(mod_fl)
} else {
  mod_real <- bayesianTrajectoryWorkflow(surv_data = bboudata::bbousurv_a %>% filter(Year > 2010), 
                                recruit_data = bboudata::bbourecruit_a %>% filter(Year > 2010),
                                niters=1,returnSamples=T)
  if(dir.exists(dirname(mod_fl))){
    saveRDS(mod_real, mod_fl)
  }
}

mod_flb <- here::here("results/test_mod_realb.rds")
if(file.exists(mod_flb)){
  mod_realb <- readRDS(mod_flb)
} else {
  disturbance <- unique(subset(bboudata::bbourecruit_a %>% filter(Year > 2010), 
                               select = Year))
  disturbance$Anthro <- 0
  disturbance$Fire_excl_anthro <- 0
  mod_realb <- bayesianTrajectoryWorkflow(
    surv_data = bboudata::bbousurv_a %>% filter(Year > 2010),
    recruit_data = bboudata::bbourecruit_a %>% filter(Year > 2010),
    disturbance = disturbance,
    niters = 10,returnSamples=T
  )
  if(dir.exists(dirname(mod_flb))){
    saveRDS(mod_realb, mod_flb)
  }
}

test_that("works with simInitial",{
  #devtools::load_all()
  # Can't compare with out disturbance in the original model
  simIni <- trajectoriesFromNational()
  
  expect_error(compareTrajectories(mod_real,
                               simInitial = simIni), "Set disturbance")
  
  
  # Works when disturbance is specified
  mod_tbl <- compareTrajectories(mod_realb, simInitial = simIni)
  
  # all components are present
  purrr::map_lgl(mod_tbl, \(x) nrow(x) > 0) %>% all() %>% 
    expect_true()
  
})

test_that("works with out simInitial", {
  
  mod_tbl <- compareTrajectories(mod_real)
  
  # sim.all is NULL when simInitial is not supplied 
  expect_null(mod_tbl$sim.all)
  
  # other components are present
  purrr::map_lgl(mod_tbl[-2], \(x) nrow(x) > 0) %>% all() %>% 
    expect_true()
  
})



test_that("exData ok in simple case with one one input scenario",{
  #exData ok in simple case with one one input scenario
  scns10 <- getScenarioDefaults(collarCount = 5, cowMult = 2, 
                                projYears = 100)
  trajs <- getSimsInitial(replicates = 2, cPars = scns10)$samples
  
  simObs8 <- simulateObservations(trajs, scns10)
  
  exDataOut <- simObs8$exData %>% 
    pivot_wider(id_cols = c("Replicate", "Year","Timestep","PopulationName"),
                names_from = "MetricTypeID",
                values_from = "Amount")
  
  # data continues to be simulated after Anthro is 100
  exDataOut %>% 
    filter(Year > 2060) %>% 
    pull(N) %>% sd() %>% {. > 1} %>% 
    expect_true()
})

test_that("sample trajectories are only returned if at least one disturbance scenario is specified", {
  expect_warning(noDist <- getSimsInitial(forceUpdate = TRUE), "a disturbance scenario must be specified")
  expect_null(noDist$samples)
})

test_that("can specify multiple disturbance scenarios", {
  scns10m <- getScenarioDefaults(collarCount = 5, cowMult = 2, 
                                 projYears = 100,iAnthro=c(0,5))
  trajs2 <- getSimsInitial(replicates = 2, cPars = scns10m)$samples
  
  # The first year will have both values for Anthro
  trajs2 %>% filter(Year == min(Year), MetricTypeID == "Anthro") %>% pull(Amount) %>% 
    range() %>% 
    expect_equal(c(0,5))
})


test_that("Warning if trajs does not include the selected disturbance scenario, and it is possible to set disturbance from trajs.", {
  scns10 <- getScenarioDefaults(collarCount = 5, cowMult = 2, 
                                projYears = 100)
  
  scns11 <- getScenarioDefaults(collarCount = 5, cowMult = 2, 
                                projYears = 100,iAnthro=7)
  trajs <- getSimsInitial(replicates = 2, cPars = scns10)$samples
  expect_warning(simulateObservations(trajs, scns11), "do not include the disturbance")
  
})

test_that("Error if trajs does not include the selected disturbance scenario, and it is possible to set disturbance from trajs.", {
  scns10m <- getScenarioDefaults(collarCount = 5, cowMult = 2, 
                                 projYears = 100,iAnthro=c(0,5))
  trajs <- getSimsInitial(replicates = 2, cPars = scns10m)$samples
  
  scns11 <- getScenarioDefaults(collarCount = 5, cowMult = 2, 
                                projYears = 100,iAnthro=7)
  expect_error(simulateObservations(trajs, scns11), "do not include the disturbance")
})

test_that("Can set disturbance from table", {
  scns10 <- getScenarioDefaults(collarCount = 5, cowMult = 2, 
                                projYears = 100)
  trajs <- getSimsInitial(data.frame(Year = seq(2009, 2017),
                                     Anthro = 5,
                                     fire_excl_anthro = 0.2),
                          replicates = 2)$samples
  
  # the traj disturbance scenario overrides scns10 disturbance scenario with a warning
  expect_warning(simObs8 <- simulateObservations(trajs, scns10), "do not include the disturbance")
  
  simObs8$exData %>% filter(MetricTypeID == "Anthro") %>% pull(Amount) %>% 
    unique() %>% expect_equal(5)
})



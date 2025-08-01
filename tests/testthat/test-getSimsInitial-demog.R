

test_that("exData ok in simple case with one one input scenario",{
  #exData ok in simple case with one one input scenario
  scns10 <- getScenarioDefaults(collarCount = 5, cowMult = 2, 
                                projYears = 100)
  simObs8 <- simulateObservations(scns10)
  
  exDataOut <- simObs8$exData %>% 
    pivot_wider(id_cols = c("Replicate", "Year","Timestep","PopulationName"),
                names_from = "MetricTypeID",
                values_from = "Amount")
  
  # data continues to be simulated after Anthro is 100, but population has collapsed
  exDataOut %>% 
    filter(Year > 2070) %>% 
    pull(N) %>% mean() %>% {. ==0} %>% 
    expect_true()

  # Visualize
  if(0){
    exDataOut %>%
      ggplot2::ggplot(ggplot2::aes(Year, N, colour = Replicate))+
      ggplot2::geom_point()+
      ggplot2::geom_point(ggplot2::aes(Year, Anthro*0.01), colour = "black")
  }
})

test_that("sample trajectories are not returned when the national model is used", {
  expect_warning(noDist <- getSimsInitial(forceUpdate = TRUE), "a disturbance scenario must be specified")
  expect_null(noDist$samples)
})

test_that("can specify multiple disturbance scenarios", {
  scns10m <- getScenarioDefaults(collarCount = 5, cowMult = 2, 
                                 projYears = 100,iAnthro=c(0,5))
  summary2 <- getSimsInitial(replicates = 2, cPars = scns10m)$summary
  
  # The first year will have both values for Anthro
  summary2 %>% filter(Year == min(Year), MetricTypeID == "N") %>% pull(AnthroID) %>% 
    range() %>% 
    expect_equal(c(0,5))
})


test_that("Warning if trajs does not include the selected disturbance scenario, and it is possible to set disturbance from trajs.", {
  scns10 <- getScenarioDefaults(collarCount = 5, cowMult = 2, 
                                projYears = 100)
  
  trajs <- getSimsInitial(replicates = 2, cPars = scns10)$samples
  expect_warning(simulateObservations(scns11,trajs), "do not include the disturbance")
  
})

test_that("Error if trajs does not include the selected disturbance scenario, and it is not possible to set disturbance from trajs.", {
  scns10m <- getScenarioDefaults(collarCount = 5, cowMult = 2, 
                                 projYears = 100,iAnthro=c(0,5))
  trajs <- getSimsInitial(replicates = 2, cPars = scns10m)$samples
  
  scns11 <- getScenarioDefaults(collarCount = 5, cowMult = 2, 
                                projYears = 100,iAnthro=7)
  expect_error(simulateObservations(scns11,trajs), "do not include the disturbance")
})

test_that("Can set disturbance from table", {
  scns10 <- getScenarioDefaults(collarCount = 5, cowMult = 2, 
                                projYears = 100)
  trajs <- getSimsInitial(data.frame(Year = seq(2009, 2017),
                                     Anthro = 5,
                                     fire_excl_anthro = 0.2),
                          replicates = 2)$samples
  
  # the traj disturbance scenario overrides scns10 disturbance scenario with a warning
  expect_warning(simObs8 <- simulateObservations(scns10,trajs), "do not include the disturbance")
  
  simObs8$exData %>% filter(MetricTypeID == "Anthro") %>% pull(Amount) %>% 
    unique() %>% expect_equal(5)
})

#TO DO: add tests trajs from bboutools models.


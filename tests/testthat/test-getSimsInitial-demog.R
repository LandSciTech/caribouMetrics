# Use saved file because this takes a long time. 
mod_fl <- here::here("results/test_mod_real.rds")
if(file.exists(mod_fl)){
  mod_real <- readRDS(mod_fl)
} else {
  mod_real <- caribouBayesianPM(surv_data = bboudata::bbousurv_a %>% filter(Year > 2010), 
                                recruit_data = bboudata::bbourecruit_a %>% filter(Year > 2010),
                                niters=1,returnSamples=T)
  if(dir.exists(dirname(mod_fl))){
    saveRDS(mod_real, mod_fl)
  }
}

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
  #devtools::load_all()
  scns10 <- getScenarioDefaults(collarCount = 5, cowMult = 2, 
                                projYears = 100)
  trajs <- subset(mod_real$result$samples,is.element(Replicate,c("x1","x2")))
  expect_warning(simulateObservations(scns11,trajs), "do not include the disturbance")
  
})

test_that("Error if trajs does not include the selected disturbance scenario, and it is not possible to set disturbance from trajs.", {
  scns10m <- getScenarioDefaults(collarCount = 5, cowMult = 2,startYear=2100)
  trajs <- subset(mod_real$result$samples,is.element(Replicate,c("x1","x2")))
  expect_error(simulateObservations(scns10m,trajs), "do not include the disturbance")
})

test_that("The trajectory can include disturbance", {
  #devtools::load_all()
  scns10 <- getScenarioDefaults(collarCount = 5, cowMult = 2, 
                                projYears = 100)
  mod_reald <- bbouMakeSummaryTable(surv_data = bboudata::bbousurv_a %>% filter(Year > 2010), 
                                 recruit_data = bboudata::bbourecruit_a %>% filter(Year > 2010),N0=1000,
                                 disturbance = data.frame(Year=seq(2010,2017),Anthro=5,fire_excl_anthro=0.2),
                                 niters=10)
  trajs <- getSimsInitial(mod_reald,replicates = 2)$samples

  # the traj disturbance scenario overrides scns10 disturbance scenario with a warning
  expect_warning(simObs8 <- simulateObservations(scns10,trajs), "do not include the disturbance")
  
  simObs8$exData %>% filter(MetricTypeID == "Anthro") %>% pull(Amount) %>% 
    unique() %>% expect_equal(5)
})


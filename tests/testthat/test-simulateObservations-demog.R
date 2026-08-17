test_that("default works", {
  scns <- getScenarioDefaults(projYears = 10, obsYears = 10, cowMult = 3,
                              collarCount = 50)
  ss <- trajectoriesFromNational(cPars=scns)
  
  expect_is(ss,"list")
  #see test-compareTrajectories-demog.R for more thorough check that outputs are as expected.
})

test_that("multiple scenarios not allowed",{
  scns <- getScenarioDefaults(data.frame(iFire = 1:2), projYears = 10, obsYears = 10)
  expect_error(simulateObservations(scns,
                                 freqStartsByYear = data.frame(Year = 2014:2023,
                                                               numStarts = 10),
                                 cowCounts = data.frame(Year = 2014:2023,
                                                        Cows = 10)),
            "cannot have multiple rows")
})

# TODO: add test for non-default popGrow table should use testPopGrowTable
# internally, do same in betaNationalPriors

test_that("collarCount and cowCount behave", {
  scns <- getScenarioDefaults(collarCount = 30, cowCount = 100, cowMult = 1)
  simObs <- simulateObservations(scns)
  
  surv_plt <- plotSurvivalSeries(simObs$simSurvObs)
  expect_is(surv_plt, "ggplot2::ggplot")
  
  # if cowCount is 100 we observe 100
  expect_true(all(simObs$simRecruitObs$Cows == 100))
  
  expect_error(getScenarioDefaults(collarCount = 30, cowCount = 100, cowMult = 2),
               "not both")
  
  # if cowMult is 2 we observe max 2*collarCount but fewer when some deaths were observed
  scns2 <- getScenarioDefaults(collarCount = 30, cowMult = 2)
  
  simObs2 <- simulateObservations(scns2)
  
  simObs2$simSurvObs %>% 
    left_join(simObs2$simRecruitObs, 
              by = join_by(PopulationName, Replicate, Year)) %>% 
    mutate(pass = Cows == (30 - MortalitiesCertain) * 2) %>% 
    pull(pass) %>% all %>% 
    expect_true()
  
  #When survival is annual, simulated survival data does not include months 1-3 of 2024. 
  #But the 2023 recruitment survey occurs in March 2024.
  expect_equal(scns2$curYear, simObs2$simSurvObs$Year %>% max())
  expect_equal(scns2$curYear+1, simObs2$simRecruitObs$Year %>% max())

  # Test with months
  # devtools::document();devtools::load_all()
  scns10 <- getScenarioDefaults(collarCount = 15, cowMult = 2)
  simObs_mon <- simulateObservations(scns10,  
                       surv_data = bboudata::bbousurv_a,
                       recruit_data = bboudata::bbourecruit_a)

  #When survival data is monthly, the survival data also includes months 1-3 of
  #2024 that will be included in the 2023 caribou year (April 2023 to March
  #2024)
  expect_equal(scns10$curYear+1, simObs_mon$simSurvObs$Year %>% max())
  expect_equal(scns10$curYear+1, simObs_mon$simRecruitObs$Year %>% max())
  
  # should be different cow counts in simulated years
  expect_gt(simObs_mon$simRecruitObs %>% filter(Year < 2016) %>% pull(Cows) %>% mean, 
            simObs_mon$simRecruitObs %>% filter(Year > 2016) %>% pull(Cows) %>% mean)
  
  # Visual:
  # simObs_mon$simRecruitObs %>% ggplot(aes(Year, Cows, colour = Replicate))+   geom_point()
  scns10 <- getScenarioDefaults(collarCount = 5, cowMult = 2, 
                                projYears = 100)
  simObs8 <- simulateObservations(scns10,  
                                  surv_data = bboudata::bbousurv_a,
                                  recruit_data = bboudata::bbourecruit_a)
  simObs8$exData %>% 
    pivot_wider(id_cols = c("Replicate", "Year","Timestep","PopulationName"),
                                 names_from = "MetricTypeID",
                                 values_from = "Amount") %>% 
  ggplot2::ggplot(ggplot2::aes(Year, N, colour = Replicate))+
    ggplot2::geom_point()+
    ggplot2::geom_point(ggplot2::aes(Year, Anthro*1), colour = "black")
  
  # if tables are supplied they should not be modified by cowCount or collarCount
  simObs3 <- simulateObservations(scns, 
                       freqStartsByYear = data.frame(Year = 2009:2023,
                                                     numStarts = 10),
                       cowCounts = data.frame(Year = 2009:2023,
                                              Cows = 10))
  
  simObs3$simRecruitObs %>% filter(Cows != 10) %>% nrow() %>% 
    {expect_true(. == 0)} 
  
  # Total n collars goes up at first and then balances with mortality
  expect_lt(
    simObs3$simSurvObs %>% filter(Year == 2010) %>% pull(StartTotal) %>% mean(),
    simObs3$simSurvObs %>% filter(Year == 2020) %>% pull(StartTotal) %>% mean()
  )
  # Visual test
  # simObs3$simSurvObs %>% 
  #   ggplot(aes(Year, StartTotal, colour = Replicate)) +
  #   geom_point()
  
  # cowMult doesn't affect cowCounts table
  simObs4 <- simulateObservations(scns2, 
                                  freqStartsByYear = data.frame(Year = 2009:2023,
                                                                numStarts = 10),
                                  cowCounts = data.frame(Year = 2009:2023,
                                                         Cows = 10))
  
  simObs4$simRecruitObs %>% filter(Cows != 10) %>% nrow() %>% 
    {expect_true(. == 0)} 
  
  # can supply just freqStartsByYear and cowMult
  simObs4b <- simulateObservations(scns2,
                                  freqStartsByYear = data.frame(Year = 2009:2023,
                                                                numStarts = 10))
  scns2$cowMult <- 10
  simObs4c <- simulateObservations(scns2, 
                                   freqStartsByYear = data.frame(Year = 2009:2023,
                                                                 numStarts = 10))
  
  # cowCounts is created from freqStartsByYear and cowMult 
  expect_true(all(
    simObs4b$simRecruitObs %>% pull(Cows) <
      simObs4c$simRecruitObs %>% pull(Cows))) 
  
  # for collarInterval 
  scns3 <- getScenarioDefaults(collarCount = 30, cowCount = 50, cowMult = NA,
                               collarInterval = 3)
  
  simObs5 <- simulateObservations(scns3)
  
  # simObs5$simSurvObs %>% ggplot(aes(Year, StartTotal, colour = Replicate)) +  geom_point()
  
  # Should go back up to 30 collars every 3 years
  simObs5$simSurvObs %>% pull(StartTotal) %>% .[seq(1,15, by = 3)] %>% 
    {. == 30} %>% all() %>% 
    expect_true()
  
  simObs5$simRecruitObs %>% filter(Cows != 50) %>% nrow() %>% 
    {expect_true(. == 0)} 
  
  # collarInterval doesn't affect tables
  simObs6 <- simulateObservations(scns3,
                                  freqStartsByYear = data.frame(Year = 2009:2023,
                                                                numStarts = 10),
                                  cowCounts = data.frame(Year = 2009:2023,
                                                         Cows = 10),
                                  collarNumYears = 20)
  
  simObs6$simSurvObs %>% group_by(PopulationName, Replicate) %>%
    mutate(recalc_starts = 10 + lag(StartTotal- MortalitiesCertain, default = 0), 
           pass = recalc_starts == StartTotal) %>% 
    pull(pass) %>% all() %>% 
    expect_true()
  
  simObs6$simRecruitObs %>% filter(Cows != 10) %>% nrow() %>% 
    {expect_true(. == 0)} 
  
  # confirm that if freqStartByYear table does skips year it still works
  simObs7 <- simulateObservations(scns,
                                  freqStartsByYear = data.frame(Year = seq(2009, 2023, by = 3),
                                                                numStarts = 10),
                                  cowCounts = data.frame(Year = 2009:2023,
                                                         Cows = 10))
  
  simObs7$simRecruitObs %>% filter(Cows != 10) %>% nrow() %>% 
    {expect_true(. == 0)} 
  
})

test_that("collarOn and Off work as expected", {
  #devtools::document();devtools::load_all()
  scns <- getScenarioDefaults(iFire = 1, projYears = 10, obsYears = 10, collarCount = 10)
  
  simObs1_12 <- simulateObservations(scns, collarOnTime = 1, collarOffTime = 12, 
                                     caribouYearStart = 4)
  #plotSurvivalSeries(simObs1_12$simSurvObs)
  expect_true(all(simObs1_12$simSurvObs$StartTotal<=scns$collarCount))
  expect_true(max(simObs1_12$simSurvObs$StartTotal)==scns$collarCount)
  expect_true(all(!is.na(simObs1_12$simSurvObs$MortalitiesCertain)))
})

#TO DO: test with trajectories from bboutools & jags models

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

test_that("exData ok in simple case with one input scenario",{
  #exData ok in simple case with one one input scenario
  scns10 <- getScenarioDefaults(collarCount = 5, cowMult = 2, 
                                projYears = 200, rSlopeMod = 2, sSlopeMod = 2)
  simObs8 <- simulateObservations(scns10)
  
  exDataOut <- simObs8$exData %>% 
    pivot_wider(id_cols = c("Replicate", "Year","Timestep","PopulationName"),
                names_from = "MetricTypeID",
                values_from = "Amount")
  
  # data continues to be simulated after Anthro is 100, but population has collapsed
  exDataOut %>% 
    filter(Year > 2200) %>% 
    pull(N) %>% mean() %>% {. == 0} %>% 
    expect_true()
  
  # Visualize
  if(0){
    exDataOut %>%
      ggplot2::ggplot(ggplot2::aes(Year, N, colour = Replicate))+
      ggplot2::geom_point()+
      ggplot2::geom_point(ggplot2::aes(Year, Anthro*0.01), colour = "black")
  }
})


test_that("Warning if trajs does not include the selected disturbance scenario, and it is possible to set disturbance from trajs.", {
  scns11 <- getScenarioDefaults(collarCount = 5, cowMult = 2, 
                                projYears = 100,iAnthro=7)
  trajs <- subset(mod_real$result$samples,is.element(Replicate,c("x1","x2")))
  expect_warning(simulateObservations(scns11,trajs), "do not include the disturbance")
  
})

test_that("Error if trajs does not include the selected disturbance scenario, and it is not possible to set disturbance from trajs.", {
  scns10m <- getScenarioDefaults(collarCount = 5, cowMult = 2, startYear=2100)
  trajs <- subset(mod_real$result$samples,is.element(Replicate,c("x1","x2")))
  expect_error(simulateObservations(scns10m,trajs), "do not include the disturbance")
})

test_that("The trajectory can include disturbance", {
  #devtools::document();devtools::load_all()
  scns10 <- getScenarioDefaults(collarCount = 5, cowMult = 2, 
                                projYears = 100)
  mod_reald <- estimateBayesianRates(surv_data = bboudata::bbousurv_a %>% filter(Year > 2010), 
                                     recruit_data = bboudata::bbourecruit_a %>% filter(Year > 2010),N0=1000,
                                     disturbance = data.frame(Year=seq(2010,2017),Anthro=5,Fire_excl_anthro=0.2),
                                     niters=15)
  
  trajs <- trajectoriesFromBayesian(mod_reald)$samples
  
  # the traj disturbance scenario overrides scns10 disturbance scenario with a warning
  expect_warning(simObs8 <- simulateObservations(scns10,trajs), "do not include the disturbance")
  
  simObs8$exData %>% filter(MetricTypeID == "Anthro") %>% pull(Amount) %>% 
    unique() %>% expect_equal(5)
})


#TO DO: test methods for combining simulated and observed data -- bboutools example
